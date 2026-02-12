/*
    integrators.c: A library for two specific numerical integrators: GSL ODEIV2 and DVERK
    Author: João Rebouças
    The purpose of integrators.c is to have a frontend for numerical algorithms.
    So far we have two options: GSL and DVERK.
 */

#include <string.h>
#include <stdbool.h>

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>


/* Wrappers for function pointer types */
typedef void (*func_dverk)(int *, double *, const double *, double *);
typedef int (*func_gsl)(double, const double *, double *, void *);


/* 
 * C interface for DVERK
 * See `dverk.f` for implementation
 */
extern void dverk_(int *n, func_dverk f,
    double *x, double *y, const double *xend, double *tol,
    int *ind, double *c, int *nw, double *w);
#define dverk dverk_

/* 
 * `dverk_opt` holds user options for DVERK
 * See the Users Guide for DVERK for more information
 * For summary:
 *   - n is the dimension of the time-dependent variable y in the ODE y' = f(y).
 *   - tol is a parameter that controls the error.
 *   - ind is an integer for setting specific configurations as described in the Users Guide.
 *     - Setting `ind == 1` represents running DVERK with default options.
 *     - Setting `ind == 2` represents running DVERK with user-defined options that are set in the array `c`.
 *     - `ind` is modified by DVERK:
 *       - If `ind == 3` after calling DVERK, the subroutine ran successfully.
 *       - If `ind<3` after calling dverk, the subroutine ran successfully.
 *   - `c` is a vector that holds additional user options. See the Users Guide for more information:
 *     - `c[0] == 1.0` means absolute error control, `c[0] == 2.0` means relative error control;
 *     - `c[2]` is the minimum step size, `c[3]` is the initial step size, `c[4]` is the characteristic time scale of the problem,  `c[5]` is the maximum step size;
 *     - `c[23]` is the number of func evals, set by DVERK;
 *   - `w` is a workspace for the subroutine.
 * NOTE: for now, c and w have arbitrary fixed capacities that I found on the examples. We might need to change the capacities
 */
#define DVERK_C_CAPACITY 24
#define DVERK_W_CAPACITY 90
typedef struct  {
    func_dverk f;
    double tol;                  
    int n;                       
    int ind;                     
    int nw;                      
    double *c;
    double *w;
} dverk_opt;

/* 
 * `gsl_opt` holds required GSL ODEIV2 structs to perform integrations.
 * See GSL manual for more information.
 */
typedef struct  {
    gsl_odeiv2_system *sys;
    gsl_odeiv2_driver *driver;
} gsl_opt;


/* Choice of integrator */
typedef enum {
    INTEGRATOR_DVERK = 0,
    INTEGRATOR_GSL   = 1,
} integrator_kind;

/*
 * `integrator_opt` is a wrapper struct that has all possible information about all integrators at once.
 * This is so that the end user doesn't need to think about which integrator to choose at allocation time.
 * It's like `integrator_opt` is a class composed by `dverk_opt` and `gsl_opt`.
 */
typedef struct {
    dverk_opt d;
    gsl_opt g;
    integrator_kind kind;
} integrator_opt;

/* 
 * `get_dverk_integrator` initializes an `integrator_opt` integrator with `.kind = INTEGRATOR_DVERK` and the `.d` field with a valid `dverk_opt`.
 */
integrator_opt get_dverk_integrator(func_dverk f, double tol, int n, int ind) {
    dverk_opt opt = { .f = f, .tol = tol, .n = n, .ind = ind, .nw = n};
    opt.c = malloc(DVERK_C_CAPACITY*sizeof(double));
    opt.w = malloc(DVERK_W_CAPACITY*sizeof(double));
    if (opt.c == NULL || opt.w == NULL) {
        fprintf(stderr, "ERROR: could not allocate memory for DVERK.");
        abort();
    }
    return (integrator_opt) { .kind = INTEGRATOR_DVERK, .d = opt };
}

/* 
 * `get_gsl_integrator` initializes an `integrator_opt` integrator with `.kind = INTEGRATOR_GSL` and the `.g` field with a valid `gsl_opt`.
 */
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
    switch (opt.kind) {
    case INTEGRATOR_GSL:
        free(opt.g.sys);
        gsl_odeiv2_driver_free(opt.g.driver);
        break;
    case INTEGRATOR_DVERK:
        free(opt.d.c);
        free(opt.d.w);
        break;
    }
}

bool integrate_dverk(double *x, double *y, const double *x_end, dverk_opt *opt) {
    dverk(&opt->n, opt->f, x, y, x_end, &opt->tol, &opt->ind, opt->c, &opt->nw, opt->w);
    return opt->ind == 3;
}

bool integrate_gsl(double *x, double *y, const double *x_end, gsl_opt *integrator_opt) {
    return gsl_odeiv2_driver_apply(integrator_opt->driver, x, *x_end, y) == GSL_SUCCESS;
}

/*
 * `integrate` performs integration of the chosen function from x to x_end, assuming initial state y.
 * x and y are going to be updated by the subroutine.
 * Returns true on success.
 */
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