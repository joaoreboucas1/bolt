#include <stdio.h>
#include "../bolt.c"

int main(void) {
    Cosmo c = {0};
    const double h = 0.67;
    const double Omega_m = 0.319;
    const double Omega_b = 0.049;
    const double A_s = 2.1e-9;
    const double n_s = 0.96;
    InitCosmo(&c, h, Omega_m, Omega_b, A_s, n_s);
    calc_background(&c);
    calc_thermo(&c);
    double z[] = {1300.0};
    Array optical_depth = get_optical_depth(z, 1);
    Array opacity = get_opacity(z, 1);
    Array visibility = get_visibility(z, 1);
    printf("At z = %g, \\kappa'(z) = %g, \\kappa(z) = %g, g(z)/H_0 = %g\n", z[0], opacity.data[0], optical_depth.data[0], visibility.data[0]);
    return 0;
}