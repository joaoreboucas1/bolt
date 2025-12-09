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

    // Arbitrary values
    double z[] = {0.5, 1.0};
    Array distances = get_luminosity_distances(z, sizeof(z)/sizeof(z[0]));
    printf("z = %.2f, D_L = %.2f Mpc || z = %.2f, D_L = %.2f Mpc\n", z[0], distances.data[0], z[1], distances.data[1]);


    return 0;
}