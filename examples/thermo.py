import bolt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

# Initialize cosmology and precompute tables
c = bolt.Cosmo(h=0.67, Omega_m=0.319, Omega_b=0.049, A_s=2.1e-9, n_s=0.96)
bolt.calc_background(c)
bolt.calc_thermo(c)

# Choose interpolation z for thermodynamical quantities
z = np.linspace(500, 2000, 1000)
optical_depth = bolt.get_optical_depth(z)
visibility = bolt.get_visibility(z)

# Make plots
fig, axs = plt.subplots(1, 2, figsize=(10, 4.5))
axs[0].plot(z, visibility, color="black")
axs[0].set_xlabel("$z$")
axs[0].set_ylabel("$g(z)$")
axs[1].plot(z, optical_depth, color="black")
axs[1].set_xlabel("$z$")
axs[1].set_ylabel("$\\kappa(z)$")
for ax in axs: 
    ax.tick_params(axis="both", which="both", direction="in", right=True, top=True)
    ax.xaxis.set_major_locator(MultipleLocator(500))
    ax.xaxis.set_minor_locator(MultipleLocator(100))
plt.savefig("Ruth_Durrer_Fig_4-1_reproduction.pdf", bbox_inches="tight")