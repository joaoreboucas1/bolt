"""
   A minimal Pythons script using Bolt.
"""

import bolt

if __name__ == "__main__":
    h = 0.67
    Omega_m = 0.319
    c = bolt.Cosmo(h, Omega_m)
    r = bolt.solve_einstein_boltzmann(c, 0.1)
    print(f"Hello, from Bolt! For h = {h:.2f}, Omega_m = {Omega_m:.3f}, \\delta_m(k=0.1, z=0) = {r.y[bolt.timesteps].delta_c:.6f}")