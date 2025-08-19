# Bolt

**Warning**: this code is made for didactical purposes only, and should not be used for serious scientific applications (yet?). For scientific applications, use [CAMB](https://github.com/cmbant/CAMB).

A library for linear cosmological perturbations. The library calculates simple transfer functions for the $\Lambda$CDM model, and can be used to calculate the matter power spectrum. We make serious simplifying assumptions which hopefully will slowly be addressed in the future.

See an example of usage in `tests.ipynb`.

`bolt` depends on [GSL](https://www.gnu.org/software/gsl/). You can download the source from the website and install it with usual commands `./configure && make && make install` to make it available system-wide. You can compile the C library using the `Makefile` we provide, possibly adjusting the GSL lib path (in my case, it was `/usr/local/lib`). Once compiled, it can be run from Python by importing the `bolt.py` library.


The core of the library is `bolt.c`, and it is interfaced with Python in `bolt.py`. The library was initially inspired in another project which is summarized in `pybolt.py`, and I wanted to reimplement it in C. Do not confuse `pybolt` with `bolt.py`. I use `pybolt` as a study to see how I want to implement things in `bolt`.

`pybolt` can be used from a Python environment with `scipy` and `numpy`.
