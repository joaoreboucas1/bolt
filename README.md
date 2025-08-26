# Bolt

**Warning**: this code is made for didactical purposes only, and should not be used for serious scientific applications (yet?). For scientific applications, use [CAMB](https://github.com/cmbant/CAMB).

A library for linear cosmological perturbations. The library calculates simple transfer functions for the $\Lambda$CDM model, and can be used to calculate the matter power spectrum. We make serious simplifying assumptions which hopefully will slowly be addressed in the future.

See an example of usage in `docs/examples.ipynb`.

## Installation

`bolt` depends on [GSL](https://www.gnu.org/software/gsl/). You can download the source from the website and install it with usual commands `./configure && make && make install` to make it available system-wide. 

To install `bolt`, first git clone this repository. Then, at its top level, use `pip` to install in editable mode:

```
    pip install -e .
```

The `setup.py` script also compiles the core C library. To compile the C library manually, we provide a `Makefile` in the `bolt/` directory alongside the source code.

You may possibly need to modify the `Makefile` to adjust the GSL lib path, since it needs to be linked with `libgsl.so` and `libgslcblas.so`. In my case, these files were located in `/usr/local/lib`.

## Code

The core of the library is `bolt.c`, which is compiled to a shared object (`bolt/build/libbolt.so`). The Python library `bolt.py` is just an interface to the C code with a few helper functions. The library was initially inspired in another [project](https://github.com/joaoreboucas1/einstein_boltzmann) which is reimplmented in `pybolt.py`. Do not confuse `pybolt` with `bolt.py`. I use `pybolt` as a study to see how I want to implement things in C for `bolt`. `pybolt` can be used from a Python environment with `scipy` and `numpy`.

## References
- [Cosmological Perturbation Equations in the Synchronous and Conformal Newtonian Gauges](https://arxiv.org/pdf/astro-ph/9506072), by C.-P. Ma and E. Bertschinger
