/*
    A minimal executable using Bolt.
*/

#include <stdio.h>
#include "../bolt.c"

int main() {
    Cosmo c = {0};
    const double h = 0.67;
    const double Omega_m = 0.319;
    InitCosmo(&c, h, Omega_m);
    Result r = solve_einstein_boltzmann(c, 0.1);
    printf("Hello, from Bolt! For h = %.2f, Omega_m = %.3f, \\delta_m(k=0.1, z=0) = %.6f\n", h, Omega_m, r.y[timesteps].delta_c);

    return 0;
}