import bolt
import numpy as np
import matplotlib.pyplot as plt

c = bolt.Cosmo(h=0.67, Omega_m=0.319, Omega_b=0.049, A_s=2.1e-9, n_s=0.96)
bolt.calc_background(c)
bolt.calc_thermo(c)
z = [1300.0]
optical_depth = bolt.get_optical_depth(z)
opacity = bolt.get_opacity(z)
visibility = bolt.get_visibility(z)
print(f"At z = {z[0]:g}, \\kappa'(z) = {opacity[0]:g}, \\kappa(z) = {optical_depth[0]:g}, g(z)/H_0 = {visibility[0]/c.H0:g}")