import bolt
import numpy as np
import matplotlib.pyplot as plt

c = bolt.Cosmo(h=0.67, Omega_m=0.319, Omega_b=0.049, A_s=2.1e-9, n_s=0.96)
bolt.calc_background(c)
z = [0.5, 1.0]
chi = bolt.get_comoving_distances(z)
print(f"z = {z[0]:.2f}, \\chi = {chi[0]:.2f} Mpc || z = {z[1]:.2f}, \\chi = {chi[1]:.2f} Mpc")