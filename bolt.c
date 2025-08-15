/*
    Bolt: a library for solving Einstein-Boltzmann system of linear cosmological perturbations
    Author: João Rebouças, August 2025
    Licence: MIT, see `LICENSE` file
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define c_in_km_s 299792.458
#define w_gamma 1.0/3.0
#define w_m 0.0
#define w_Lambda -1.0
const int timesteps = 3999;

typedef struct {
    // Input parameters
    double h;
    double Omega_m;
    // Derived parameters
    double Omega_Lambda;
    double H0;
    double rho_crit;
    double Omega_r;
    double a_eq;
} Cosmo;

Cosmo InitCosmo(double h, double Omega_m) {
    const double H0 = 100*h/c_in_km_s;
    const double Omega_r = 2.5e-5;
    return (Cosmo) {
        .h = h,
        .Omega_m = Omega_m,
        .H0 = H0, // 1/Mpc
        .rho_crit = 3.0*H0*H0, // 1/Mpc^2
        .Omega_r = Omega_r,
        .Omega_Lambda = 1.0 - Omega_m - Omega_r,
        .a_eq = Omega_r/Omega_m
    };
}

// Background densities
double rho_m(Cosmo c, double a) {
    return c.rho_crit*c.Omega_m*pow(a, -3.0);
}

double rho_gamma(Cosmo c, double a) {
    return c.rho_crit*c.Omega_r*pow(a, -4.0);
}

double rho_lambda(Cosmo c, double a) {
    (void) a;
    return c.rho_crit*c.Omega_Lambda;
}

double rho_tot(Cosmo c, double a) {
    return c.rho_crit*(c.Omega_r*pow(a, -4.0) + c.Omega_m*pow(a, -3.0) + c.Omega_Lambda);
}

double P(Cosmo c, double a) {
    return c.rho_crit*(w_gamma*c.Omega_r*pow(a, -4.0) + w_m*c.Omega_m*pow(a, -3.0) + w_Lambda*c.Omega_Lambda);
}

double H_curly(Cosmo c, double a) {
    return a*sqrt(rho_tot(c, a)/3.0);
}

// Einstein-Boltzmann system
double scale_factor_horizon_entry(Cosmo c, double k) {
    (void) c;
    (void) k;
    // Find the value of `a` such that k = H_curly
    // Python implementation:
    // return np.exp(root_scalar(lambda loga: H_curly(np.exp(loga)) - k, x0=-15).root)
    printf("scale_factor_horizon_entry: not implemented");
    exit(1);
    return 0.0;
}

typedef struct {
    double delta_c;
    double theta_c;
    double delta_gamma;
    double theta_gamma;
    double Phi;
} Perturbations;

void Perturbations_Scale(Perturbations *p, double s) {
    p->delta_c *= s;
    p->theta_c *= s;
    p->delta_gamma *= s;
    p->theta_gamma *= s;
    p->Phi *= s;
}

Perturbations Perturbations_Add(Perturbations a, Perturbations b) {
    return (Perturbations) {
        .delta_c = a.delta_c + b.delta_c,
        .theta_c = a.theta_c + b.theta_c,
        .delta_gamma = a.delta_gamma + b.delta_gamma,
        .theta_gamma = a.theta_gamma + b.theta_gamma,
        .Phi = a.Phi + b.Phi,
    };
}

Perturbations dy_da(Cosmo c, Perturbations y, double a, double k) {
    const double k2 = pow(k, 2.0);
    const double H = H_curly(c, a);
    const double rho_m_now = rho_m(c, a);
    const double rho_gamma_now = rho_gamma(c, a);
    const double sigma_gamma = 0.0; // TODO: implement anisotropic stress, for now it's zero
    
    const double Phi_prime   = -H*y.Phi + 0.5*pow(a,2.0)*(rho_m_now*y.theta_c + rho_gamma_now*(1 + w_gamma)*y.theta_gamma)/k2;
    const double delta_c_prime = -y.theta_c + 3.0*Phi_prime;
    const double theta_c_prime = -H*y.theta_c + k2 * y.Phi;
    const double delta_gamma_prime = -4.0*y.theta_gamma/3.0 + 4.0*Phi_prime;
    const double theta_gamma_prime = k2*(y.delta_gamma/4.0 - sigma_gamma) + k2*y.Phi;

    Perturbations p = (Perturbations) {
        .delta_c = delta_c_prime,
        .theta_c = theta_c_prime,
        .delta_gamma = delta_gamma_prime,
        .theta_gamma = theta_gamma_prime,
        .Phi = Phi_prime
    };
    Perturbations_Scale(&p, 1.0/(a*H));
    return p;
}

Perturbations dy_dloga(Cosmo c, Perturbations y, double loga, double k) {
    const double a = exp(loga);
    Perturbations y_prime = dy_da(c, y, a, k);
    Perturbations_Scale(&y_prime, a);
    return y_prime;
}

// odeint(dy_dloga, y, loga_int, k);
void odeint(Cosmo c, Perturbations f(Cosmo c, Perturbations y, double loga, double k), Perturbations *y, double *loga_int, double k) {
    Perturbations y_now, y_prime;
    const double dloga_int = loga_int[2] - loga_int[1];
    for (int i = 0; i < timesteps; i++) {
        y_now = y[i];
        y_prime = f(c, y_now, loga_int[i], k);
        Perturbations_Scale(&y_prime, dloga_int);
        y[i+1] = Perturbations_Add(y_now, y_prime);
    }
}

typedef struct {
    Perturbations *y;
    double *a;
    int timesteps;
} Result;

Result integrate(Cosmo c, double k) {
    // Initial conditions
    const double tau_ini = 3e-4;
    const double a_ini = c.H0*sqrt(c.Omega_r)*tau_ini + pow(c.H0,2.0)*c.Omega_m*pow(tau_ini,2.0)/4.0;
    const double C = 1.0; // Arbitrary scale of adiabatic initial conditions
    const double Phi_ini = 4*C/3; // Primordial curvature perturbation
    const double delta_gamma_ini = -2*Phi_ini;
    const double theta_gamma_ini = 0.5*pow(k,2.0)*tau_ini*Phi_ini;
    const double delta_c_ini = 3*delta_gamma_ini/4;
    const double theta_c_ini = theta_gamma_ini;

    Perturbations y_ini = (Perturbations) {
        .Phi = Phi_ini,
        .delta_gamma = delta_gamma_ini,
        .theta_gamma = theta_gamma_ini,
        .delta_c = delta_c_ini,
        .theta_c = theta_c_ini,
    };
    
    const int n_data = timesteps + 1;
    Perturbations *y = malloc(n_data*sizeof(Perturbations));
    double *a_int = malloc(n_data*sizeof(double));
    if (y == NULL || a_int == NULL) {
        printf("ERROR: could not allocate memory for integration\n");
        exit(1);
    }

    y[0] = y_ini;

    // Defining scale factor grid for integration
    double loga_int[n_data];
    const double dloga_int = -log(a_ini)/timesteps;
    for (int i = 0; i < n_data; i++) {
        loga_int[i] = log(a_ini) + i*dloga_int;
        a_int[i] = exp(loga_int[i]);
    }
    odeint(c, dy_dloga, y, loga_int, k);

    return (Result) {.y = y, .a = a_int, .timesteps=timesteps};
}