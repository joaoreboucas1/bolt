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
    calc_transfers(&c);
    // printf("For k = %e, z = %f, \\delta_c = %e", ks[num_k-1], bg.z[timesteps], transfer_functions[num_k-1][timesteps].delta_c);
    double z[] = {0.0, 0.5, 1.0};
    double k[] = {1e-3, 1e-2, 1e-1};
    Array tk = get_matter_tk(k, 3, z, 3);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            printf("For k = %.3f, z = %.3f, \\delta_c = %.3f\n", k[i], z[j], tk.data[j*3 + i]);
        }
    }

    return 0;
}