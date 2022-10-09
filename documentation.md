project: SciTools 
src_dir: ./src
output_dir: ./doc
exclude: lebedev_weights.F90 
	units_inc.f90
exclude_dir: ./src/math/lebedev/
author: Michael Schueler
author_description: The SciTools library: A collection of mathematical and numerical tools for scientific applications.
project_github: https://github.com/michaelschueler/scitools
project_download: https://github.com/michaelschueler/scitools/releases/latest
github: https://github.com/michaelschueler


Description
-----------------

**SciTools** is a collection of mathematical and numerical tools for scientific applications.


It includes tools for:

* interpolation ([[scitools_bsplines]]), finite difference ([[scitools_fornberb]]) schemes for differentiation, Gregory quadrature ([[scitools_gregory]]), Runge-Kutta integration ([[scitools_rungekutta]]).

* Gauss-Legendre quadrature ([[scitools_gausslegendre]]), Filon quadrature ([[scitools_filon]]), Lebedev quadrature for integrals over the unit sphere ([[scitools_lebedev_quad]]).

* simplified calls to BLAS/LAPACK routines ([[scitools_linalg]]).

* defining of constants, tools for in/output ([[scitools_utils]]), a simple timer ([[scitools_time]]).

* reading and writing data from / to text and binary files ([[scitools_binutils]]).

* computing matrix exponentials ([[scitools_matrixexp]]) and time evolution of quantum mechanical systems ([[scitools_evol]]).

Optional features:

* Simple parallel MPI distribution scheme for one-dimensional arrays ([[scitools_array1d_dist]]).

* Interfaces for reading/writing data from/to HDF5 files ([[scitools_hdf5_utils]]).

The shared library `libscitools` can be used to compile other programs.

Official Releases
-----------------

The **current stable release** can be [downloaded
on GitHub](https://github.com/michaelschueler/scitools/releases/latest). The
documentation for the current version is included in the the release.

A list of all past releases, links to their documentation, and the
change log can be found on the
[releases page](https://github.com/michaelschueler/scitools/releases/index.html).

Installation
------------

### Dependencies ###
* [cmake](https://cmake.org), version >= 2.8
* A version of BLAS/LAPACK 
* [hdf5](https://www.hdfgroup.org/solutions/hdf5/) - optional
* MPI - optional
* The code documentation tool [FORD](https://github.com/cmacmackin/ford) is
used for the source code documentation.

### Compilation ###

We use the [cmake](https://cmake.org) build system to generate Makefiles. We recommend creating a configuration script 
with all options:

```
# From the root directory .../scitools/
# Create configure script
vim configure.sh
```

The basic configure script has the following structure:

```
FC=[Fortran compiler] \
cmake \
	-DCMAKE_INSTALL_PREFIX = install_path \
	-DCMAKE_BUILD_TYPE=[Release|Debug] \
	-DCMAKE_INCLUDE_PATH=/opt/local/include \
	-DCMAKE_LIBRARY_PATH=/opt/local/lib \
	-DCMAKE_Fortran_FLAGS="[optimization flags|debug flags]" \
	..
```

In the first line, the Fortran compiler is defined. The install directory (for instance `/home/opt`) is defined by the CMake variable `CMAKE_INSTALL_PREFIX`. Debugging tools are switched on by setting `CMAKE_BUILD_TYPE=Debug`; otherwise, all assertions and sanity checks are turned off. 

As the next step create a build directory

```
# From the root directory .../scitools/
# Create build directory
mkdir build
cd build
```

and run the configure script:

```
sh ../configure.sh
```

After successful configuration (which generates the Makefiles), compile the library by

```
make
```

and install via

```
make install
```
#### BLAS/LAPACK ####

A version of BLAS / LAPACK is required. In addition to specify the location of BLAS/LAPACK in the `CMAKE_LIBRARY_PATH`, the user should specify their version of BLAS/LAPACK via the `lapackblas_libraries` variable in the configure script. This is required in most cases (exception: compiler wrappers that automatically link against BLAS/LAPACK). 

For example, if you are using Intel MKL, add

	cmake -Dlapackblas_libraries="-mkl"

or similar (refer to the MKL link advisor). For the OpenBLAS library, specify

	cmake -Dlapackblas_libraries="-lopenblas"

Other implementations of BLAS/LAPACK work similarly.

#### HDF5 ####

To enable hdf5 support use the `hdf5` option in cmake

    cmake -Dhdf5=ON

Make sure that the HDF5 library, including the Fortran interface library, is installed on your systems and the library path is included in `CMAKE_LIBRARY_PATH`. Also add the include path where `hdf5.mod` is located to `CMAKE_INCLUDE_PATH`.

#### MPI ####

To enable MPI support use the `mpi` option in cmake

    cmake -Dmpi=ON

You will need to specify an MPI Fortran compiler, for example

    FC=mpif90 \
    cmake \ 
    	...	

