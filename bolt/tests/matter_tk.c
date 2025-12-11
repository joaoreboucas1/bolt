/*
    A minimal executable using Bolt.
*/

#include <stdio.h>
#include "../bolt.c"

int main() {
    Cosmo c = {0};
    const double h = 0.67;
    const double Omega_m = 0.319;
    const double Omega_b = 0.049;
    const double A_s = 2.1e-9;
    const double n_s = 0.96;

    InitCosmo(&c, h, Omega_m, Omega_b, A_s, n_s);
    calc_background(&c);
    solve_einstein_boltzmann(c, 0.1);
    printf("Hello, from Bolt! For h = %.2f, Omega_m = %.3f, \\delta_m(k=0.1, z=0) = %.6f\n", h, Omega_m, result[timesteps].delta_c);

    return 0;
}