"""
    Bolt.py: Python bindings for the Bolt library
    Author: João Rebouças, August 2025
    Licence: MIT, see `LICENSE` file
"""

import ctypes
from typing import *
import numpy as np

libbolt = ctypes.CDLL("./libbolt.so")

def InitCosmo(h, Omega_m):
    return libbolt.InitCosmo(h, Omega_m)

class Cosmo(ctypes.Structure):
    _fields_ = [
        ("h", ctypes.c_double),
        ("Omega_m", ctypes.c_double),
        ("Omega_Lambda", ctypes.c_double),
        ("H0", ctypes.c_double),
        ("rho_crit", ctypes.c_double),
        ("Omega_r", ctypes.c_double),
        ("a_eq", ctypes.c_double),
    ]

    def __init__(self, h: float, Omega_m: float) -> Self:
        self = libbolt.InitCosmo(h, Omega_m)

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}: \n" + f"  - h = {self.h:.4f}\n" + f"  - Omega_m = {self.Omega_m:.4f}"
    
    def rho_m(self, a):
        return libbolt.rho_m(self, a)
    
    def rho_gamma(self, a: float) -> float:
        return libbolt.rho_gamma(self, a)

    def rho_lambda(self, a: float) -> float:
        return libbolt.rho_lambda(self, a)

    def rho_tot(self, a: float) -> float:
        return libbolt.rho_tot(self, a)

    def P(self, a: float) -> float:
        return libbolt.P(self, a)

    def H_curly(self, a: float) -> float:
        return libbolt.H_curly(self, a)

    def scale_factor_horizon_entry(self, k: float) -> float:
        return libbolt.scale_factor_horizon_entry(self, k)

class Perturbations(ctypes.Structure):
    _fields_ = [
        ("delta_c", ctypes.c_double),
        ("theta_c", ctypes.c_double),
        ("delta_gamma", ctypes.c_double),
        ("theta_gamma", ctypes.c_double),
        ("Phi", ctypes.c_double),
    ]

    def as_np_array(self):
        return np.array([getattr(self, field[0]) for field in self._fields_])
    
    def from_list(l) -> Self:
        return Perturbations(*l)
    
    def __repr__(self):
        return f"{self.__class__.__name__}: \n" + "\n".join([f"  - {field[0]}: {getattr(self, field[0])}" for field in self._fields_])

class Result(ctypes.Structure):
    _fields_ = [
        ("y", ctypes.POINTER(Perturbations)),
        ("a", ctypes.POINTER(ctypes.c_double)),
        ("timesteps", ctypes.c_int),
    ]

    def as_arrays(self):
        a = np.ctypeslib.as_array(self.a, shape=(self.timesteps+1,))
        y = np.ctypeslib.as_array(self.y, shape=(self.timesteps+1,))
        return a, y

# int dy_dloga(double loga, Perturbations y, Perturbations *y_prime, void *params) {

def solve_einstein_boltzmann(c: Cosmo, k: float) -> Result:
    return libbolt.solve_einstein_boltzmann(c, k)

def integrate(c: Cosmo, k: float) -> Result:
    # Just for backwards compatibility
    return solve_einstein_boltzmann(c, k)


libbolt.rho_m.argtypes = (Cosmo, ctypes.c_double)
libbolt.rho_m.restype = ctypes.c_double
libbolt.rho_gamma.argtypes = (Cosmo, ctypes.c_double)
libbolt.rho_gamma.restype = ctypes.c_double
libbolt.rho_lambda.argtypes = (Cosmo, ctypes.c_double)
libbolt.rho_lambda.restype = ctypes.c_double
libbolt.rho_tot.argtypes = (Cosmo, ctypes.c_double)
libbolt.rho_tot.restype = ctypes.c_double
libbolt.P.argtypes = (Cosmo, ctypes.c_double)
libbolt.P.restype = ctypes.c_double
libbolt.H_curly.argtypes = (Cosmo, ctypes.c_double)
libbolt.H_curly.restype = ctypes.c_double
libbolt.scale_factor_horizon_entry.argtypes = (Cosmo, ctypes.c_double)
libbolt.scale_factor_horizon_entry.restype = ctypes.c_double
libbolt.InitCosmo.argtypes = [ctypes.c_double, ctypes.c_double]
libbolt.InitCosmo.restype = Cosmo
libbolt.dy_dloga.argtypes = [Cosmo, Perturbations, ctypes.c_double, ctypes.c_double]
libbolt.dy_dloga.restype = Perturbations
libbolt.solve_einstein_boltzmann.argtypes = [Cosmo, ctypes.c_double]
libbolt.solve_einstein_boltzmann.restype = Result

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    c = libbolt.InitCosmo(0.7, 0.3)
    k = 0.01
    r = libbolt.solve_einstein_boltzmann(c, 0.01)
    delta_c = [r.y[i].delta_c for i in range(r.timesteps)]
    a = np.array([r.a[i] for i in range(r.timesteps)])
    plt.loglog(a, np.abs(delta_c))
    plt.show()
    # data = np.array(r.y)
    