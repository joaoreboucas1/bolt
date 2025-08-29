# Bolt

> [!Warning]
> This code is made for didactical purposes only, and should not be used for serious scientific applications (yet?). For scientific applications, use [CAMB](https://github.com/cmbant/CAMB).

A library for linear cosmological perturbations. The library calculates simple transfer functions for the $\Lambda$CDM model, and can be used to calculate the CDM power spectrum. We make serious simplifying assumptions which hopefully will be slowly addressed in the future.

See an example of usage in `docs/examples.ipynb`.

## Installation

> [!Tip]
> I recommend reading the installation steps before running any command.

> [!Important]
> If any installation step fails, please open an issue describing your problem in detail, providing error messages.

1. `bolt` depends on [GSL](https://www.gnu.org/software/gsl/). If you already have a GSL installation, you can skip to the next step. The following command downloads the GSL source code and performs a system-wide installation:

```console
    wget https://mirror.ibcp.fr/pub/gnu/gsl/gsl-latest.tar.gz
    tar -xvf gsl-latest.tar.gz # This command creates a folder named `gsl-X.Y`, where X.Y is the GSL version. In the following, I will use version 2.8, but it may be different for you
    cd gsl-2.8/
    ./configure
    make
    make install
```

> [!Note]
> You may want to install GSL in another folder if, for instance, you don't have permission to install GSL system-wide (e.g. in a cluster). For this purpose, change the `./configure` command to `./configure --prefix=/some/path/to/gsl/`, where `/some/path/to/gsl` is the path you want to install it. In this case, the GSL headers will be located in `/some/path/to/gsl/include/` and the compiled librbaries will be located in `/some/path/to/gsl/lib/`.

> [!Note]
> Clusters may already provide GSL installations. You can check the output of the `module avail` command. If you see GSL installations in the output, you can run the command `module load gsl-X.Y`. Double-check where are the headers and compiled libraries; if you cannot find them, contact the system administrator. If your cluster does not provide a GSL installation, you can ask the system administrator to install it system-wide.

If compilation and installation went successfully, you should find the GSL headers in `/usr/local/include/gsl/` and the compiled libraries (`libgsl.so.28` and `libgslcblas.so`) in `/usr/local/lib/`.

2. To install `bolt`, first git clone this repository. Then, at its top level, use `pip` to install in editable mode:

```
    pip install -e .
```

The `setup.py` script compiles the core C library. To compile the C library manually, we provide a `Makefile` in the `bolt/` directory alongside the source code. You may possibly need to modify the `Makefile` to adjust the GSL lib path according to Step 1. In particular, modify the variables `GSL_INC_PATH` and `GSL_LIB_PATH` according to your needs.

To test the installation, you can run the `docs/examples.ipynb` Jupyter notebook.

## Code

The core of the library is `bolt.c`, containing all of the background functions and linear perturbation equations. The C code is compiled to a shared object (`bolt/build/libbolt.so`). The Python library `bolt.py` is an interface to the C code with a few helper functions.

The library was initially inspired in another [project](https://github.com/joaoreboucas1/einstein_boltzmann) which is reimplmented in `pybolt.py`. Do not confuse `pybolt` with `bolt.py`. I use `pybolt` as a study to see how I want to implement things in C for `bolt`. `pybolt` can be used from a Python environment with `scipy` and `numpy`.

## References
- [Cosmological Perturbation Equations in the Synchronous and Conformal Newtonian Gauges](https://arxiv.org/pdf/astro-ph/9506072), by C.-P. Ma and E. Bertschinger
