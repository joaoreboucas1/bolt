/*
    Bolt: a library for solving Einstein-Boltzmann system of linear cosmological perturbations
    Author: João Rebouças, August 2025
    Licence: MIT, see `LICENSE` file
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>

#define c_in_km_s 299792.458
#define omega_gammah2 2.473e-5 // From PDG 2020 astrophysical constants
#define w_gamma 1.0/3.0
#define w_m 0.0
#define w_Lambda -1.0
#define L_MAX_NU 100
const int timesteps = 999;

// TODO: maybe put a_int, As and ns here
typedef struct {
    // Input parameters
    double h;
    double Omega_m;
    // Derived parameters
    double Omega_Lambda;
    double H0;
    double rho_crit;
    double Omega_gamma;
    double a_eq;
} Cosmo;

void InitCosmo(Cosmo *c, double h, double Omega_m) {
    const double H0 = 100*h/c_in_km_s;
    const double Omega_gamma = omega_gammah2/h/h;
    *c = (Cosmo) {
        .h = h,
        .Omega_m = Omega_m,
        .H0 = H0, // 1/Mpc
        .rho_crit = 3.0*H0*H0, // 1/Mpc^2
        .Omega_gamma = Omega_gamma,
        .Omega_Lambda = 1.0 - Omega_m - Omega_gamma,
        .a_eq = Omega_gamma/Omega_m
    };
}

// Background densities
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

typedef struct {
    double delta_c;
    double theta_c;
    double delta_gamma;
    double theta_gamma;
    double Phi;
    double tau; // TODO: build a background table
} Perturbations;

static_assert(sizeof(Perturbations) == 6*sizeof(double), "Exhaustive handling of Perturbations in Perturbations_Scale");
void Perturbations_Scale(Perturbations *p, double s) {
    p->delta_c *= s;
    p->theta_c *= s;
    p->delta_gamma *= s;
    p->theta_gamma *= s;
    p->Phi *= s;
    p->tau *= s;
}

static_assert(sizeof(Perturbations) == 6*sizeof(double), "Exhaustive handling of Perturbations in Perturbations_Add");
Perturbations Perturbations_Add(Perturbations a, Perturbations b) {
    Perturbations p = (Perturbations) {
        .delta_c = a.delta_c + b.delta_c,
        .theta_c = a.theta_c + b.theta_c,
        .delta_gamma = a.delta_gamma + b.delta_gamma,
        .theta_gamma = a.theta_gamma + b.theta_gamma,
        .Phi = a.Phi + b.Phi,
        .tau = a.tau + b.tau,
    };
    return p;
}

typedef struct {
    Cosmo c;
    double k;
} ODE_Params;

// int deriv(double t, const double y[], double dydt[], void *Params)
static_assert(sizeof(Perturbations) == 6*sizeof(double), "Exhaustive handling of Perturbations in dy_da");
int dy_da(double a, Perturbations y, Perturbations *y_prime, void *params) {
    ODE_Params ode_params = *(ODE_Params*) params;
    const Cosmo c = ode_params.c;
    const double k = ode_params.k;
    const double k2 = k*k;
    const double a2 = a*a;
    const double H = H_curly(c, a);
    const double rho_m_now = rho_m(c, a);
    const double rho_gamma_now = rho_gamma(c, a);
    const double sigma_gamma = 0.0; // TODO: implement anisotropic stress, for now it's zero
    
    const double Phi_prime   = -H*y.Phi + 0.5*a2*(rho_m_now*y.theta_c + rho_gamma_now*(1.0 + w_gamma)*y.theta_gamma)/k2;
    const double delta_c_prime = -y.theta_c + 3.0*Phi_prime;
    const double theta_c_prime = -H*y.theta_c + k2 * y.Phi;
    const double delta_gamma_prime = -4.0*y.theta_gamma/3.0 + 4.0*Phi_prime;
    const double theta_gamma_prime = k2*(y.delta_gamma/4.0 - sigma_gamma) + k2*y.Phi;
    const double tau_prime = 1.0;

    *y_prime = (Perturbations) {
       .delta_c = delta_c_prime,
       .theta_c = theta_c_prime,
       .delta_gamma = delta_gamma_prime,
       .theta_gamma = theta_gamma_prime,
       .Phi = Phi_prime,
       .tau = tau_prime
   };

    Perturbations_Scale(y_prime, 1.0/(a*H));
    return GSL_SUCCESS;
}

int dy_dloga(double loga, Perturbations y, Perturbations *y_prime, void *params) {
    const double a = exp(loga);
    dy_da(a, y, y_prime, params);
    Perturbations_Scale(y_prime, a);
    return GSL_SUCCESS;
}

int dy_dloga_gsl(double loga, const double y[], double y_prime[], void *params) {
    return dy_dloga(loga, *(Perturbations*)y, (Perturbations*)y_prime, params);
}

typedef struct {
    Perturbations *y;
    double *loga;
    int timesteps;
} Result;

// odeint(dy_dloga, y, loga_int, k);
Result odeint(Cosmo c, Perturbations y_ini, double loga_ini, double k) {
    const double dloga_int = -loga_ini/timesteps;
    ODE_Params params = (ODE_Params) {
        .c = c,
        .k = k
    };

    gsl_odeiv2_system sys = {
        .function = dy_dloga_gsl,
        .jacobian = NULL,
        .dimension = (int)(sizeof(Perturbations)/sizeof(double)),
        .params = (void*)&params
    };
    const double hstart = dloga_int;
    const double epsabs = 1e-3;
    const double epsrel = 0.0;
    gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, hstart, epsabs, epsrel);

    const int n_data = timesteps + 1;
    // TODO: I'm returning newly allocated memory to the caller. Is this a good idea?
    Perturbations *y = malloc(n_data*sizeof(Perturbations));
    double *loga_int = malloc(n_data*sizeof(double));
    if (loga_int == NULL) {
        printf("ERROR: could not allocate memory for integration\n");
        exit(1);
    }
    Perturbations state = y_ini;
    y[0] = state;

    double loga = loga_ini;
    loga_int[0] = loga_ini;
    for (int i = 0; i < timesteps; i++) {
        const double loga_next = loga + dloga_int;
        gsl_odeiv2_driver_apply(driver, &loga, loga_next, (double*)&state);
        // TODO: maybe consider fixed step?
        // gsl_odeiv2_driver_apply_fixed_step(driver, &loga, dloga_int, 1, (double*)&state);
        y[i + 1] = state;
        loga_int[i + 1] = loga;
    }
    gsl_odeiv2_driver_free(driver);
    return (Result) {.loga = loga_int, .y = y, .timesteps = timesteps};
}

static_assert(sizeof(Perturbations) == 6*sizeof(double), "Exhaustive handling of Perturbations in initial_conditions");
Perturbations initial_conditions(Cosmo c, double k, double tau_ini) {
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
        .tau = tau_ini
    };
    return y_ini;
}

Result solve_einstein_boltzmann(Cosmo c, double k) {
    const double tau_ini = 3e-4;
    const double a_ini = c.H0*sqrt(c.Omega_gamma)*tau_ini + c.H0*c.H0*c.Omega_m*tau_ini*tau_ini/4.0;
    Perturbations y_ini = initial_conditions(c, k, tau_ini);
    return odeint(c, y_ini, log(a_ini), k);
}