# SciTools

**SciTools** is a Fortran library of mathematical and numerical tools for
scientific applications. It builds a single library, `libscitools`, whose
modules (all prefixed `scitools_`) can be used from other Fortran programs.

SciTools is self-contained and serves as the numerical foundation of
higher-level simulation codes.

## Features

| Area | Modules |
|------|---------|
| Basics: kinds, constants, units, I/O helpers, assertions | `def`, `utils`, `debug` |
| Linear algebra: simplified BLAS/LAPACK interfaces | `linalg` |
| Iterative linear solvers | `cg_solver`, `bicgstab`, `gmresm`, `idrs` |
| Quadrature | `gausslegendre`, `gauss_lobatto`, `gregory`, `filon`, `quadpack`, `quadrature`, `lebedev_quad` (unit sphere) |
| Interpolation and differentiation | `bsplines`, `vector_bsplines`, `fornberg` (finite differences) |
| ODEs and time evolution | `rungekutta`, `matrixexp`, `evol` |
| Special functions and polynomials | `special`, `legendre`, `hermite` |
| Root finding, sorting, random numbers | `root`, `sorting`, `random` |
| Timing | `time` |
| File I/O | `binutils` (text/binary), `hdf5_utils` (HDF5, optional) |
| Parallelization | `array1d_dist`: distribution of 1D arrays over MPI ranks (optional) |
| Physics helpers | `laserpulse`, `gausspulse` |

## Installation

### Requirements

* [CMake](https://cmake.org) >= 3.10
* A Fortran compiler (Fortran 2008) and a C compiler
* A BLAS/LAPACK implementation (e.g. OpenBLAS, Intel MKL, Apple Accelerate)
* Optional: [HDF5](https://www.hdfgroup.org/solutions/hdf5/) with the Fortran interface
* Optional: an MPI implementation with Fortran bindings
* Optional: [FORD](https://github.com/Fortran-FOSS-Programmers/ford) to generate the API documentation

### Building

SciTools uses an out-of-source CMake build. A template configure script is
provided in `configure.sh.example`:

```sh
# from the scitools root directory
cp configure.sh.example configure.sh   # edit to match your system
mkdir build && cd build
sh ../configure.sh
make
make install
```

A minimal configure script looks like this:

```sh
FC=gfortran \
cmake \
    -DCMAKE_INSTALL_PREFIX=$HOME/opt \
    -DCMAKE_BUILD_TYPE=Release \
    -Dlapackblas_libraries="-lopenblas" \
    -DCMAKE_Fortran_FLAGS="-cpp -O3 -fPIC" \
    ..
```

Setting `CMAKE_BUILD_TYPE=Debug` enables assertions and sanity checks, which
are disabled in release builds. If dependencies are installed in
non-standard locations, add them to `CMAKE_INCLUDE_PATH` and
`CMAKE_LIBRARY_PATH`.

### Build options

| Option | Default | Description |
|--------|---------|-------------|
| `-Dlapackblas_libraries="..."` | unset | Link flags for BLAS/LAPACK, e.g. `-lopenblas`, `-mkl`, `-framework Accelerate`. Required unless your compiler wrapper links BLAS/LAPACK automatically. |
| `-Dhdf5=ON` | `OFF` | Build `scitools_hdf5_utils`. HDF5 is located with CMake's `find_package(HDF5 COMPONENTS Fortran)`. |
| `-Dmpi=ON` | `OFF` | Build `scitools_array1d_dist` and enable MPI code paths. MPI is located with `find_package(MPI COMPONENTS Fortran)`. |
| `-Dmodernquad=ON` | `OFF` | Use the modern quadrature backend in `scitools_quadrature`. |
| `-Denable_static=ON` | `OFF` | Build a static library instead of a shared one. |

### Installed files

`make install` places

* the library in `<prefix>/lib/libscitools.{so,dylib,a}`,
* the Fortran module files in `<prefix>/include/scitools/scitools_*.mod`.

## Using SciTools in other projects

Add the module directory to the include path and link against `libscitools`
(and the same BLAS/LAPACK, HDF5 and MPI libraries used to build it):

```sh
gfortran -I$HOME/opt/include/scitools myprog.f90 -L$HOME/opt/lib -lscitools -lopenblas
```

```fortran
program example
  use scitools_def, only: dp
  use scitools_linalg, only: eye, inv
  implicit none
  real(dp) :: A(3,3), Ainv(3,3)

  A = 2.0_dp * eye(3)
  Ainv = inv(A)
  print *, Ainv(1,1)   ! 0.5
end program example
```

In CMake projects, SciTools can either be used from an existing installation
(pass the install prefix and link `libscitools`) or be built as part of the
project with `FetchContent`:

```cmake
include(FetchContent)
FetchContent_Declare(scitools
    GIT_REPOSITORY https://github.com/michaelschueler/scitools.git
    GIT_TAG        <release tag>)
set(lapackblas_libraries "-lopenblas" CACHE STRING "" FORCE)  # plus hdf5/mpi as needed
FetchContent_MakeAvailable(scitools)
target_link_libraries(myprog PRIVATE scitools)
target_include_directories(myprog PRIVATE ${scitools_BINARY_DIR}/include)
```

The SciTools build options (`lapackblas_libraries`, `hdf5`, `mpi`,
`enable_static`) are ordinary cache variables and can be set before
`FetchContent_MakeAvailable`.

## Documentation

API documentation is generated from the source code with FORD (version 7 or
later):

```sh
ford documentation.md
```

The HTML documentation is written to `./doc`; open `doc/index.html`.

## Third-party code and acknowledgements

Parts of SciTools are based on, or adapted from, the following works. We
thank the original authors. Details and the full license texts are given in
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) and the [licenses](licenses/)
directory.

| SciTools module | Original work | Authors | License |
|-----------------|---------------|---------|---------|
| `linalg`, `utils`, `sorting`, `special`, `debug` | [fortran-utils](https://github.com/certik/fortran-utils); parts based on code by J. E. Pask (LLNL) | O. Čertík | MIT |
| `special` | `SPHJ`, `SPHY`, `MSTA1`, `MSTA2` from `specfun.f` (via SciPy) | S. Zhang, J. Jin, *Computation of Special Functions* (Wiley, 1996) | Permission granted to use in programs |
| `bsplines` | [bspline-fortran](https://github.com/jacobwilliams/bspline-fortran), based on NIST CMLIB and SLATEC | J. Williams; NIST; C. de Boor | BSD-3-Clause / public domain |
| `gausslegendre` | [quadrature-fortran](https://github.com/jacobwilliams/quadrature-fortran), including SLATEC `DGAUS8` | J. Williams; SLATEC | BSD-3-Clause / public domain |
| `quadpack` | QUADPACK | R. Piessens, E. de Doncker-Kapenger, C. Überhuber, D. Kahaner | Public domain |
| `hermite`, `filon` | Routines by J. Burkardt (`hermite` partly after S. Elhay and J. Kautsky) | J. Burkardt | MIT |
| `idrs` | IDR(s) solver | P. Sonneveld, M. B. van Gijzen, SIAM J. Sci. Comput. 31, 1035 (2008) | <!-- TODO --> to be confirmed |
| `gmresm` | GMRES(m) with hookstep for Newton-Krylov methods | A. P. Willis <!-- TODO: confirm origin --> | to be confirmed |
| `lebedev_weights` | Lebedev quadrature grids | V. I. Lebedev, D. N. Laikov, Doklady Mathematics 59, 477 (1999) <!-- TODO: data source --> | |
| `gauss_lobatto` | Gauss-Lobatto nodes and weights | <!-- TODO: data source --> | |

Algorithms implemented from the literature include the matrix exponential by
scaling and squaring (N. J. Higham, SIAM J. Matrix Anal. Appl. 26, 1179 (2005)),
BiCGSTAB (H. A. van der Vorst, SIAM J. Sci. Stat. Comput. 13, 631 (1992)), the
Fornberg finite-difference weights (B. Fornberg, Math. Comp. 51, 699 (1988)) and
Gregory quadrature.

## License

SciTools is licensed under the Mozilla Public License 2.0; see
[LICENSE.md](LICENSE.md). Third-party code keeps its original license; see
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
