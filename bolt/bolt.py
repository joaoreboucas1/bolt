"""
    Bolt.py: Python bindings for the Bolt library
    Author: João Rebouças, August 2025
    Licence: MIT, see `LICENSE` file
"""

import os
import ctypes
from typing import *
import numpy as np

libc = ctypes.CDLL("libc.so.6")
libc.free.argtypes = (ctypes.c_void_p,)
libc.free.restype = None

path = "/".join(__file__.split("/")[:-1])
if not os.path.exists(f"{path}/build/libbolt.so"):
    raise ImportError(f"Bolt tried to look for `libbolt.so` at `{path}/libbolt.so`, but it does not exist.")

try:
    libbolt = ctypes.CDLL(f"{path}/build/libbolt.so")
except OSError as e:
    raise ImportError(f"Bolt tried to load `libbolt.so`, but we got an error: \n{e}\n. Make sure compilation went fine; otherwise, you can create an issue in Github.")

# TODO: the .so does not export macros. How can we get constants?
with open(f"{path}/bolt.c", "r") as f:
    for line in f.read().splitlines():
        if line.startswith("#define timesteps"):
            timesteps = int(line.split()[-1])

class Cosmo(ctypes.Structure):
    _fields_ =  [
        ("h", ctypes.c_double),
        ("Omega_m", ctypes.c_double),
        ("Omega_b", ctypes.c_double),
        ("A_s", ctypes.c_double),
        ("n_s", ctypes.c_double),
        ("Omega_Lambda", ctypes.c_double),
        ("Omega_c", ctypes.c_double),
        ("H0", ctypes.c_double),
        ("rho_crit", ctypes.c_double),
        ("Omega_gamma", ctypes.c_double),
        ("a_eq", ctypes.c_double),
    ]

    def __init__(self, h: float, Omega_m: float, Omega_b: float, A_s: float, n_s: float) -> Self:
        libbolt.InitCosmo(ctypes.byref(self), h, Omega_m, Omega_b, A_s, n_s)

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}:\n" + "\n".join([f"  - {field[0]}: {getattr(self, field[0])}" for field in self._fields_[:5]])
    
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

class Array(ctypes.Structure):
    _fields_ = [("data", ctypes.POINTER(ctypes.c_double)),
                ("len", ctypes.c_int)]


def calc_background(c: Cosmo):
    libbolt.calc_background(ctypes.byref(c))

def get_comoving_distances(z_values):
    z_array = np.asarray(z_values, dtype=np.float64)
    result = libbolt.get_comoving_distances(z_array.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), len(z_array))
    # Convert C array to numpy array
    d_L = np.ctypeslib.as_array(result.data, shape=(result.len,)).copy()
    
    # Free C-allocated memory
    libc.free(result.data)
    
    return d_L

def solve_einstein_boltzmann(c: Cosmo, k: float):
    # TODO: libbolt.solve_einstein_boltzmann returns a pointer to a global variable, is this a good idea?
    return libbolt.solve_einstein_boltzmann(c, k)

def integrate(c: Cosmo, k: float):
    # Backwards compatibility
    return solve_einstein_boltzmann(c, k)

# C function bindings
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
libbolt.InitCosmo.argtypes = [ctypes.POINTER(Cosmo), ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double]
libbolt.calc_background.argtypes = [ctypes.POINTER(Cosmo)]
libbolt.calc_background.restype = None
libbolt.get_comoving_distances.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_size_t]
libbolt.get_comoving_distances.restype = Array
# TODO: interface dy_da
# libbolt.dy_dloga.argtypes = [Cosmo, Perturbations, ctypes.c_double, ctypes.c_double]
# libbolt.dy_dloga.restype = Perturbations
libbolt.solve_einstein_boltzmann.argtypes = [Cosmo, ctypes.c_double]
libbolt.solve_einstein_boltzmann.restype = ctypes.POINTER(Perturbations)    