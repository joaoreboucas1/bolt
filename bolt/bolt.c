/*
 *  Bolt: a library for solving Einstein-Boltzmann system of linear cosmological perturbations
 *  Author: João Rebouças, August 2025
 *  Licence: MIT, see `LICENSE` file
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_errno.h>

#include "integrators.c"

// -------------------- Constants --------------------

#define c_in_km_s 299792.458

/*   
 *  Species:
 *      - `gamma` refers to photons;
 *      - `m` refers to matter (CDM + baryons);
 *      - `Lambda` refers to the cosmological constant \Lambda;
 *      - `w_<species>` refers to the equation of state of species <species>;
 *      - `omega_<species>` refers to the current density parameter $ \Omega_i = \rho_i(z=0)/\rho_tot(z=0) $ of species <species>;
 */

#define omega_gammah2 2.473e-5 // From PDG 2020 astrophysical constants
#define w_gamma 1.0/3.0
#define w_m 0.0
#define w_Lambda -1.0

/*
 *  Time-stepping:
 *      - Perturbation equations are integrated over $ln(a)$ and require background predictions
 *      - We choose an initial conformal time $\tau_\mathrm{ini} = 3x10^(-4) Mpc$
 *      - We calculate the initial value of "a" using a radiation-dominated approximation (requires cosmological parameters)
 *      - We choose linearly separated steps on `loga` until `loga = 0`
 *      - We precompute a table of background quantities (a, tau, rho_i, H) calculated at the `loga` values
 *      - These tables are used to generate interpolators
 */

 // NOTE: the number of time points, including the initial and final times, is `timesteps + 1` 
#define timesteps 999
#define num_loga (timesteps+1)
#define num_a num_loga
const double tau_ini = 3e-4;

// -------------------- Model Interface --------------------

/*
 *  `Cosmo` is a struct that holds input parameters as well as useful derived parameters that can be reused.
 */

typedef struct {
    // Input parameters
    double h;
    double Omega_m;
    double Omega_b;
    double A_s;
    double n_s;

    // Derived parameters set by `InitCosmo`
    double Omega_Lambda;
    double Omega_c;
    double H0;
    double rho_crit;
    double Omega_gamma;
    double a_eq;
} Cosmo;

// Global variable that integrator can access
static Cosmo c_global;

void InitCosmo(Cosmo *c, double h, double Omega_m, double Omega_b, double A_s, double n_s) {
    const double H0 = 100*h/c_in_km_s;
    const double Omega_gamma = omega_gammah2/h/h;
    *c = (Cosmo) {
        // Input parameters
        .h = h,
        .Omega_m = Omega_m,
        .Omega_b = Omega_b,
        .A_s = A_s,
        .n_s = n_s,

        // Derived parameters
        .Omega_Lambda = 1.0 - Omega_m - Omega_gamma,
        .Omega_c = Omega_m - Omega_b,
        .H0 = H0,              // 1/Mpc
        .rho_crit = 3.0*H0*H0, // 1/Mpc^2
        .Omega_gamma = Omega_gamma,
        .a_eq = Omega_gamma/Omega_m
    };

    // NOTE: we have a global variable that holds the current cosmo object
    c_global = *c;
}

// -------------------- Background functions --------------------

/*
 *  Quantities:
 *      - `rho` is the energy density in units of 1/Mpc^2 (because 8*pi*G=1, see "Constants" section);
 *      - `P` is the pressure in the same units (c = 1);
 *      - `H_curly` is the conformal Hubble factor.
 */

double rho_m(Cosmo c, double a) {
    return c.rho_crit*c.Omega_m/a/a/a;
}

double rho_gamma(Cosmo c, double a) {
    return c.rho_crit*c.Omega_gamma/a/a/a/a;
}

double rho_lambda(Cosmo c, double a) {
    (void) a;
    return c.rho_crit*c.Omega_Lambda;
}

double rho_tot(Cosmo c, double a) {
    return c.rho_crit*(c.Omega_gamma/a/a/a/a + c.Omega_m/a/a/a + c.Omega_Lambda);
}

double P(Cosmo c, double a) {
    return c.rho_crit*(w_gamma*c.Omega_gamma/a/a/a/a + w_m*c.Omega_m/a/a/a + w_Lambda*c.Omega_Lambda);
}

double H_curly(Cosmo c, double a) {
    return a*sqrt(rho_tot(c, a)/3.0);
}

/*
 *  `Background` is a type for the library's background tables.
 *  The integration of perturbations requires the computation of several background functions.
 *  The ultimate goal of the library is to calculate perturbation equations for several values of the Fourier wavenumber $k$.
 *  We don't need to recalculate the background for each Fourier wavenumber, so we just do it once and save the results in background tables.
 *  Then, we define GSL interpolators that can be reused in the perturbation equations as well as in the outputs of the library.
 */

typedef struct {
    double loga[timesteps+1];
    double a[timesteps+1];
    double z[timesteps+1];
    double tau[timesteps+1];
    double H[timesteps+1];
} Background;

Background bg;
gsl_interp *tau_interpolator = NULL;
gsl_interp *H_interpolator = NULL;

double inverse_NormHubble_gsl_integration(double z, void* params) {
    // NOTE: NormHubble = H(z)/H_0
    Cosmo c = *(Cosmo*)params;
    double a = 1.0/(1.0 + z);
    double H = H_curly(c, a)/a;
    return c.H0/H;
}

// Computes the background table and initializes GSL interpolators
void calc_background(Cosmo *c) {
    const double a_ini = c->H0*sqrt(c->Omega_gamma)*tau_ini + c->H0*c->H0*c->Omega_m*tau_ini*tau_ini/4.0;
    bg.a[0] = a_ini;
    bg.loga[0] = log(a_ini);
    bg.z[0] = 1.0/a_ini - 1.0;
    bg.tau[0] = tau_ini;
    bg.H[0] = H_curly(*c, a_ini)/a_ini;
    const double dloga_int = -bg.loga[0]/timesteps;
    
    // NOTE: \tau(z) - \tau(z_ini) = \int_{z}^{z_ini} dz/H(z)
    double integral, abserr;
    double cumulative_integral = 0.0;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(timesteps);
    gsl_function integrand;
    integrand.function = &inverse_NormHubble_gsl_integration;
    integrand.params = (void*)c;
    for (size_t i = 1; i < timesteps+1; ++i) {
        bg.loga[i]  = bg.loga[i-1] + dloga_int;
        bg.a[i]     = exp(bg.loga[i]);
        bg.z[i]     = 1.0/bg.a[i] - 1.0;
        bg.H[i]     = H_curly(*c, bg.a[i])/bg.a[i];

        // Performing \tau integration
        int status = gsl_integration_qag(
            &integrand,
            bg.z[i], 
            bg.z[i-1], 
            0.0,  // abserr
            1e-5, // relerr
            timesteps, 
            GSL_INTEG_GAUSS21, 
            w,  
            &integral, 
            &abserr);
        if (status != GSL_SUCCESS) {
            fprintf(stderr, "ERROR: could not integrate background tau");
            exit(1);
        }
        cumulative_integral += integral;
        bg.tau[i] = tau_ini + 1.0/c->H0*cumulative_integral; // 
    }

    if (tau_interpolator == NULL) tau_interpolator = gsl_interp_alloc(gsl_interp_cspline, timesteps+1);
    if (H_interpolator == NULL) H_interpolator = gsl_interp_alloc(gsl_interp_cspline, timesteps+1);
    int status = gsl_interp_init(tau_interpolator, bg.loga, bg.tau, timesteps+1);
    if (status != GSL_SUCCESS) {
        fprintf(stderr, "ERROR: could not initialize tau interpolator\n");
        exit(1);
    }
    status = gsl_interp_init(H_interpolator, bg.loga, bg.H, timesteps+1);
    if (status != GSL_SUCCESS) {
        fprintf(stderr, "ERROR: could not initialize Hubble factor interpolator\n");
        exit(1);
    }
    gsl_integration_workspace_free(w);
}

// Interface with NumPy `np.array`
typedef struct {
    double *data;
    size_t len;
} Array;

// API function that gets the comoving distance in Mpc to a sources with redshift `z`, an array with length `z_len`.
// Must be called after calling `calc_background()`. 
Array get_comoving_distances(double *z, size_t z_len) {
    double *data = malloc(z_len*sizeof(double));
    if (data == NULL) {
        fprintf(stderr, "ERROR: could not allocate memory for comoving distances\n");
        exit(1);
    }
    gsl_interp_accel *accel = gsl_interp_accel_alloc();
    for (size_t i = 0; i < z_len; ++i) {
        float a = 1.0/(1.0+z[i]);
        // TODO: error when out of bounds
        float result = gsl_interp_eval(tau_interpolator, bg.loga, bg.tau, log(a), accel); // tau(a)
        data[i] = bg.tau[timesteps] - result;
    }
    gsl_interp_accel_free(accel);
    return (Array) { .data = data, .len = z_len };
}

// -------------------- Perturbations --------------------

/*
 *  The `Perturbations` struct fields:
 *      - `delta` denotes the density contrast;
 *      - `theta` denotes the velocity divergence, $ \theta = k*v $;
 *      - `Phi` is the Newtonian gauge gravitational potential;
 *  The suffix of each field is the species, e.g. `delta_c` is the CDM density contrast.
 */

typedef struct {
    double delta_c;
    double theta_c;
    double delta_gamma;
    double theta_gamma;
    double Phi;
} Perturbations;

// Convenience function to multiply `Perturbation` by some factor
static_assert(sizeof(Perturbations) == 5*sizeof(double), "Exhaustive handling of Perturbations in Perturbations_Scale");
void Perturbations_Scale(Perturbations *p, double s) {
    p->delta_c *= s;
    p->theta_c *= s;
    p->delta_gamma *= s;
    p->theta_gamma *= s;
    p->Phi *= s;
}

// Convenience function to add two `Perturbation`s
static_assert(sizeof(Perturbations) == 5*sizeof(double), "Exhaustive handling of Perturbations in Perturbations_Add");
Perturbations Perturbations_Add(Perturbations a, Perturbations b) {
    Perturbations p = (Perturbations) {
        .delta_c = a.delta_c + b.delta_c,
        .theta_c = a.theta_c + b.theta_c,
        .delta_gamma = a.delta_gamma + b.delta_gamma,
        .theta_gamma = a.theta_gamma + b.theta_gamma,
        .Phi = a.Phi + b.Phi,
    };
    return p;
}

// TODO: make `k_global` _Thread_local
double k_global;

// Core function for the derivatives of scalar perturbations.
// Functions for interfacing with specific integrator (e.g. GSL, DVERK) can be implemented in terms of this function
static_assert(sizeof(Perturbations) == 5*sizeof(double), "Exhaustive handling of Perturbations in dy_da");
int dy_da(double a, Perturbations y, Perturbations *y_prime) {    
    const double k2 = k_global*k_global;
    const double a2 = a*a;

    // TODO: we could use the background interpolators for these calculations
    const double H = H_curly(c_global, a);
    const double rho_m_now = rho_m(c_global, a);
    const double rho_gamma_now = rho_gamma(c_global, a);
    const double sigma_gamma = 0.0; // TODO: implement anisotropic stress, for now it's zero
    
    const double Phi_prime   = -H*y.Phi + 0.5*a2*(rho_m_now*y.theta_c + rho_gamma_now*(1.0 + w_gamma)*y.theta_gamma)/k2;
    const double delta_c_prime = -y.theta_c + 3.0*Phi_prime;
    const double theta_c_prime = -H*y.theta_c + k2 * y.Phi;
    const double delta_gamma_prime = -4.0*y.theta_gamma/3.0 + 4.0*Phi_prime;
    const double theta_gamma_prime = k2*(y.delta_gamma/4.0 - sigma_gamma) + k2*y.Phi;

    *y_prime = (Perturbations) {
       .delta_c = delta_c_prime,
       .theta_c = theta_c_prime,
       .delta_gamma = delta_gamma_prime,
       .theta_gamma = theta_gamma_prime,
       .Phi = Phi_prime,
    };
    Perturbations_Scale(y_prime, 1.0/(a*H));
    return GSL_SUCCESS;
}

// Helper function that converts dy_da into dy_dloga
int dy_dloga(double loga, Perturbations y, Perturbations *y_prime) {
    const double a = exp(loga); // TODO: can we use the interpolation table?
    dy_da(a, y, y_prime);
    Perturbations_Scale(y_prime, a);
    return GSL_SUCCESS;
}

// Helper function to convert `dy_dloga` to the GSL format
int dy_dloga_gsl(double loga, const double y[], double y_prime[], void *params) {
    (void)params;
    return dy_dloga(loga, *(Perturbations*)y, (Perturbations*)y_prime);
}

// Helper function to convert `dy_dloga` to the DVERK format
void dy_dloga_dverk(int *ndim, double *loga, const double *y, double *y_prime) {
    (void)ndim;
    dy_dloga(*loga, *(Perturbations*)y, (Perturbations*)y_prime);
}

static_assert(sizeof(Perturbations) == 5*sizeof(double), "Exhaustive handling of Perturbations in initial_conditions");
Perturbations initial_conditions(Cosmo c, double k) {
    (void) c;
    const double k2 = k*k;
    const double C = 1.0; // Arbitrary scale of adiabatic initial conditions
    const double Phi_ini = 4.0*C/3.0; // Primordial curvature perturbation
    const double delta_gamma_ini = -2.0*Phi_ini;
    const double theta_gamma_ini = 0.5*k2*tau_ini*Phi_ini;
    const double delta_c_ini = 3.0*delta_gamma_ini/4.0;
    const double theta_c_ini = theta_gamma_ini;

    Perturbations y_ini = (Perturbations) {
        .Phi = Phi_ini,
        .delta_gamma = delta_gamma_ini,
        .theta_gamma = theta_gamma_ini,
        .delta_c = delta_c_ini,
        .theta_c = theta_c_ini,
    };
    return y_ini;
}

// Solves the Einstein-Boltzmann system for a single value of `k`.
// Stores the result in `result`, which must be an array of `Perturbations` of size `timesteps+1`
void solve_einstein_boltzmann(Cosmo cosmo, double k, Perturbations *result) {
    int n_dim = (int)(sizeof(Perturbations)/sizeof(double));
    double tol = 1e-3;
    
    integrator_opt opt1 = get_gsl_integrator(dy_dloga_gsl, tol, n_dim);
    // TODO: allocate `c` and `w` on the heap
    double c[DVERK_C_CAPACITY];
    double w[DVERK_W_CAPACITY];
    integrator_opt opt2 = get_dverk_integrator(dy_dloga_dverk, tol, n_dim, 1, c, w);
    
    const Perturbations y_ini = initial_conditions(cosmo, k);
    Perturbations state = y_ini;
    result[0] = state;
    
    k_global = k; // NOTE: the derivative function only knows about k_global so we must set it
    double loga = bg.loga[0];
    const double dloga_int = -bg.loga[0]/timesteps;
    for (int i = 0; i < timesteps; i++) {
        const double loga_next = loga + dloga_int;
        integrate(&loga, (double*)&state, &loga_next, &opt2);
        result[i + 1] = state;
    }
    
    integrator_free(opt1);
    integrator_free(opt2);
}

// Bolt has a default array of `k` values in h/Mpc
#define logk_min -3.0
#define logk_max 0
#define num_logk 10
#define num_k num_logk
#define dlogk (logk_max - logk_min)/(num_logk-1)
static_assert(num_logk > 1);

// Buffer to store the result of `compute_transfers`
Perturbations transfer_functions[num_logk][timesteps+1];
double ks[num_logk];

void calc_transfers(Cosmo *c) {
    // TODO: this can be parallelized
    for (size_t i = 0; i < num_logk; ++i) {
        double logk = logk_min + i*dlogk;
        double k = pow(10.0, logk);
        ks[i] = k;
        solve_einstein_boltzmann(*c, k, transfer_functions[i]);
    }
}

Array get_matter_tk(double *k, size_t k_len, double *z, size_t z_len) {
    gsl_interp2d *matter_tk_interp = gsl_interp2d_alloc(gsl_interp2d_bilinear, num_k, num_loga);
    double matter_tk[num_k*num_loga];
    
    // TODO: use gsl_interp2d_idx and gsl_interp2d_set
    // TODO: the `matter_tk_interp` object could be computed during `calc_transfers` like the background
    for (size_t i = 0; i < num_k; ++i) {
        for (size_t j = 0; j < num_loga; ++j) {
            matter_tk[j*num_k + i] = transfer_functions[i][j].delta_c;
        }
    }
    
    int status = gsl_interp2d_init(
        matter_tk_interp,
        ks,
        bg.loga,
        matter_tk,
        num_k,
        num_loga
    );

    if (status != GSL_SUCCESS) {
        fprintf(stderr, "ERROR: could not initialize gsl_interp2d\n");
        exit(1);
    }
    
    double *data = calloc(k_len * z_len, sizeof(double));
    if (data == NULL) {
        fprintf(stderr, "ERROR: could not allocate memory for interpolation of matter tk\n");
        exit(1);
    }

    gsl_interp_accel *xaccel = gsl_interp_accel_alloc();
    gsl_interp_accel *yaccel = gsl_interp_accel_alloc();
    
    for (size_t i = 0; i < k_len; ++i) {
        for (size_t j = 0; j < z_len; ++j) {
            double loga = log(1.0/(1.0+z[j]));
            int status = gsl_interp2d_eval_extrap_e(matter_tk_interp, ks, bg.loga, matter_tk, k[i], loga, xaccel, yaccel, &data[j*k_len + i]);
            if (status != GSL_SUCCESS) printf("WARNING: in matter tk interpolation, gsl_interp_2d gave an error");
        }
    }

    gsl_interp2d_free(matter_tk_interp);
    gsl_interp_accel_free(xaccel);
    gsl_interp_accel_free(yaccel);

    return (Array) { .data = data, .len = k_len*z_len };
}