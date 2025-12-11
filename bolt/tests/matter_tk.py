"""
   A minimal Python script using Bolt.
"""

import bolt

if __name__ == "__main__":
    h = 0.67
    Omega_m = 0.319
    Omega_b = 0.049
    A_s = 2.1e-9
    n_s = 0.96

    c = bolt.Cosmo(h, Omega_m, Omega_b, A_s, n_s)
    bolt.calc_background(c)
    result = bolt.solve_einstein_boltzmann(c, 0.1)
    print(f"Hello, from Bolt! For h = {h:.2f}, Omega_m = {Omega_m:.3f}, \\delta_m(k=0.1, z=0) = {result[bolt.timesteps].delta_c:.6f}")