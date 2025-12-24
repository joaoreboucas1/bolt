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
                ("len", ctypes.c_size_t)]

def calc_background(c: Cosmo):
    libbolt.calc_background(ctypes.byref(c))

def get_comoving_distances(z_values):
    z_array = np.asarray(z_values, dtype=np.float64)
    result = libbolt.get_comoving_distances(z_array.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), len(z_array))
    # Convert C array to numpy array
    d_L = np.ctypeslib.as_array(result.data, shape=(result.len,)).copy()
    
    # Free C-allocated memory
    libc.free(ctypes.cast(result.data, ctypes.c_void_p))
    
    return d_L

def calc_transfers(c: Cosmo):
    libbolt.calc_transfers(ctypes.byref(c))

def get_matter_tk(k_values, z_values):
    z_array = np.asarray(z_values, dtype=np.float64)
    k_array = np.asarray(k_values, dtype=np.float64)
    
    result = libbolt.get_matter_tk(k_array.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), len(k_array), z_array.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), len(z_array))
    
    # Convert C buffer to a 1D numpy array, reshape into (z_len, k_len) as stored in C, then transpose
    flat = np.ctypeslib.as_array(result.data, shape=(result.len,)).copy()
    tk = flat.reshape((len(z_array), len(k_array))).T

    # Free C-allocated memory
    libc.free(result.data)

    return tk

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
libbolt.calc_transfers.argtypes = [ctypes.POINTER(Cosmo)]
libbolt.calc_transfers.restype = None
libbolt.get_comoving_distances.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_size_t]
libbolt.get_comoving_distances.restype = Array
libbolt.get_matter_tk.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_size_t, ctypes.POINTER(ctypes.c_double), ctypes.c_size_t]
libbolt.get_matter_tk.restype = Array