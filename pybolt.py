"""
    pyBolt: initial explorations for the Boly library
    Author: João Rebouças, August 2025
    Licence: MIT, see `LICENSE` file
"""

import numpy as np
from scipy.integrate import odeint
from scipy.optimize import root_scalar

# Constants
c_in_km_s = 299_792.458
w_gamma = 1/3
w_m = 0
w_Lambda = -1
Omega_r = 2.5e-5

# Free and derived parameters
h = 0.7
H0 = 100*h/c_in_km_s # 1/Mpc
rho_crit = 3*H0**2 # 1/Mpc^2
Omega_m = 0.3
Omega_Lambda = 1 - Omega_m - Omega_r
a_eq = Omega_r/Omega_m

class Cosmo:
    def __init__(self, h: float, Omega_m: float):
        self.h = h
        self.Omega_m = Omega_m
        self.Omega_r = Omega_r
        self.H0 = 100*h/c_in_km_s
        self.rho_crit = 3*H0**2
        self.Omega_Lambda = 1 - Omega_m - Omega_r
        self.a_eq = Omega_r/Omega_m

    def __repr__(self):
        return "Cosmology: \n" + f"  - h = {self.h:.4f}\n" + f"  - Omega_m = {self.Omega_m:.4f}"

    # Densities (8*\pi*G*\rho)
    def rho_m(self, a: float):
        return self.rho_crit*self.Omega_m*a**-3

    def rho_gamma(self, a: float):
        return self.rho_crit*self.Omega_r*a**-4

    def rho_lambda(self, a: float):
        return self.rho_crit*self.Omega_Lambda

    def rho_tot(self, a: float):
        return self.rho_crit*(self.Omega_r*a**-4 + self.Omega_m*a**-3 + self.Omega_Lambda)

    def P(self, a: float):
        return self.rho_crit*(w_gamma*self.Omega_r*a**-4 + w_m*self.Omega_m*a**-3 + w_Lambda*self.Omega_Lambda)

    def H_curly(self, a: float):
        return a*np.sqrt(self.rho_tot(a)/3)

    def scale_factor_horizon_entry(self, k: float):
        # Find the value of `a` such that k = H_curly
        return np.exp(root_scalar(lambda loga: self.H_curly(np.exp(loga)) - k, x0=-15).root)

# Derivative of the system with respect to scale factor
def dy_da(y: np.ndarray, a: float, c: Cosmo, k: float) -> np.ndarray:
    delta_m, theta_m, delta_gamma, theta_gamma, Phi = y
    k2 = k**2
    H = c.H_curly(a)
    rho_m_now = c.rho_m(a)
    rho_gamma_now = c.rho_gamma(a)
    sigma_gamma = 0 # TODO: implement anisotropic stress, for now it's zero
    Phi_prime   = -H*Phi + 0.5*a**2*(rho_m_now*theta_m + rho_gamma_now*(1 + w_gamma)*theta_gamma)/k2
    delta_m_prime = -theta_m + 3*Phi_prime
    theta_m_prime = -H*theta_m + k2*Phi
    delta_gamma_prime = -4*theta_gamma/3 + 4*Phi_prime
    theta_gamma_prime = k2*(delta_gamma/4 - sigma_gamma) + k2*Phi
    return np.array([delta_m_prime, theta_m_prime, delta_gamma_prime, theta_gamma_prime, Phi_prime])/(a*H)

def dy_dloga(y: np.ndarray, loga: float, c: Cosmo, k: float) -> np.ndarray:
    a = np.exp(loga)
    deriv = dy_da(y, a, c, k)
    return a*deriv

def integrate(c: Cosmo, k: float) -> tuple[np.ndarray, np.ndarray]:
    # Initial conditions
    tau_ini = 3e-4
    a_ini = c.H0*np.sqrt(c.Omega_r)*tau_ini + c.H0**2*c.Omega_m*tau_ini**2/4
    C = 1 # Arbitrary scale of adiabatic initial conditions
    Phi_ini = 4/3*C # Primordial curvature perturbation
    delta_gamma_ini = -2*Phi_ini
    theta_gamma_ini = 0.5*k**2*tau_ini*Phi_ini
    delta_c_ini = 3*delta_gamma_ini/4
    theta_c_ini = theta_gamma_ini
    y0 = [delta_c_ini, theta_c_ini, delta_gamma_ini, theta_gamma_ini, Phi_ini]

    # Defining scale factor grid for integration
    loga_int = np.linspace(np.log(a_ini), np.log(1), 4000)
    a = np.exp(loga_int)
    result = odeint(dy_dloga, y0, t=loga_int, args=(c, k))
    return a, result

def solve_system_for_ks(c: Cosmo, ks: np.typing.ArrayLike) -> tuple[np.ndarray, list[np.ndarray]]:
    results = [integrate(c, k)[1] for k in ks]
    for k in ks:
        a, result = integrate(c, k)
        results.append(result)
    return a, results