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
    bolt.calc_transfers(c)

    zz = [0, 0.5, 1.0]
    kk = [1e-3, 1e-2, 1e-1]

    tk = bolt.get_matter_tk(kk, zz)
    for i, k in enumerate(kk):
        for j, z in enumerate(zz):
            print(f"For k = {k:.3f}, z = {z:.3f}, \\delta_z = {tk[i][j]:.3f}")
