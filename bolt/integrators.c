/*
    integrators.c: A library for two specific numerical integrators
    Author: João Rebouças
    The purpose of integrators.c is to have a frontend for numerical algorithms.
    So far we have two options: GSL and DVERK.
 */

#include <string.h>
#include <stdbool.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>


/* Wrappers for function pointer types */
typedef void (*func_dverk)(int *, double *, double *, double *);
typedef int (*func_gsl)(double, const double *, double *, void *);

/* DVERK interface */
extern void dverk_(int *n, func_dverk fcn,
    double *x, double *y, const double *xend, double *tol,
    int *ind, double *c, int *nw, double *w);
#define dverk dverk_

typedef enum {
    INTEGRATOR_DVERK = 0,
    INTEGRATOR_GSL   = 1,
} integrator_kind;


/* dverk_opt holds user options for DVERK, see the Users Guide for DVERK for more information */
/* n is the dimension of the ODE */
/* tol is an error controlling parameter */
/* The standard DVERK options are chosen by passing ind = 1 */
/* ind is modified by dverk. If ind==3 after calling dverk, the subroutine ran successfully. */
/* c is a vector that holds additional user options, see the Users Guide for DVERK for more information */
/* w is a workspace for the subroutine */
/* NOTE: for now, c and w have arbitrary fixed capacities that I found on the examples. */
#define DVERK_C_CAPACITY 24
#define DVERK_W_CAPACITY 90
typedef struct  {
    func_dverk f;
    double tol;                  // Abs error
    int n;                       // Dimension of the system
    int ind;                     // `ind = 1` means no options to be used; After the `dverk` subroutine, `ind < 0` indicate errors. `ind = 2` indicates the use of `c`
    int nw;                      // Usually equal to `n`
    double c[DVERK_C_CAPACITY];  // c[0] = 1.0 means absolute error control, c[0] = 2.0 means relative error control;
                                 // c[2] = hmin, c[3] = hstart, c[4] is the characteristic time scale,  c[5] = hmax, c[23] is the number of func evals
    double w[DVERK_W_CAPACITY];
} dverk_opt;

/* gsl_opt holds ODEsystem definitions for GSL. */
/* See GSL manual for more information. */
typedef struct  {
    gsl_odeiv2_system *sys;
    gsl_odeiv2_driver *driver;
} gsl_opt;

/* integrator_opt is a wrapper that has all possible information about all integrators at once */
typedef struct {
    dverk_opt d;
    gsl_opt g;
    integrator_kind kind;
} integrator_opt;

/* Initializes a DVERK integrator */
integrator_opt get_dverk_integrator(func_dverk f, double tol, int n, int ind, double c[DVERK_C_CAPACITY], double w[DVERK_W_CAPACITY]) {
    dverk_opt opt = { .f = f, .tol = tol, .n = n, .ind = ind, .nw = n};
    memcpy(&opt.c, &c[0], DVERK_C_CAPACITY*sizeof(double));
    memcpy(&opt.w, &w[0], DVERK_W_CAPACITY*sizeof(double));
    return (integrator_opt) { .kind = INTEGRATOR_DVERK, .d = opt };
}

/* Initializes a GSL integrator */
integrator_opt get_gsl_integrator(func_gsl f, double tol, int n) {
    double hstart = 0.1;
    double epsrel = 0.0;
    integrator_opt opt = {0};
    opt.kind = INTEGRATOR_GSL;
    opt.g.sys = malloc(sizeof(gsl_odeiv2_system));
    *(opt.g.sys) = (gsl_odeiv2_system) { .function = f, .jacobian = NULL, .dimension = n, .params = NULL };
    opt.g.driver = gsl_odeiv2_driver_alloc_y_new(opt.g.sys, gsl_odeiv2_step_rkf45, hstart, tol, epsrel);
    return opt;
}

/* Frees the integrator */
void integrator_free(integrator_opt opt) {
    if (opt.kind == INTEGRATOR_GSL) {
        free(opt.g.sys);
        gsl_odeiv2_driver_free(opt.g.driver);
    }
}

bool integrate_dverk(double *x, double *y, const double *x_end, dverk_opt *opt) {
    dverk(&opt->n, opt->f, x, y, x_end, &opt->tol, &opt->ind, opt->c, &opt->nw, opt->w);
    return opt->ind == 3;
}

bool integrate_gsl(double *x, double *y, const double *x_end, gsl_opt *integrator_opt) {
    return gsl_odeiv2_driver_apply(integrator_opt->driver, x, *x_end, y) == GSL_SUCCESS;
}

/* Performs integration of the chosen function from x to x_end, assuming initial state y */
/* x and y are going to be updated by the subroutine */
/* Returns true on success */
bool integrate(double *x, double *y, const double *x_end, integrator_opt *opt) {
    switch (opt->kind) {
    case INTEGRATOR_DVERK:
        return integrate_dverk(x, y, x_end, &opt->d);
        break;
    case INTEGRATOR_GSL:
        return integrate_gsl(x, y, x_end, &opt->g);
        break;
    default:
        fprintf(stderr, "ERROR: unknown integrator_kind in integrate");
        abort();
        break;
    }
}