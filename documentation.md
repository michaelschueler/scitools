project: SciTools
src_dir: ./src
output_dir: ./doc
exclude: **/*_inc.f90
exclude_dir: ./src/math/lebedev
             ./src/math/lobatto
docmark_alt: ~
author: Michael Schueler
author_description: The SciTools library: A collection of mathematical and numerical tools for scientific applications.
project_github: https://github.com/michaelschueler/scitools
project_download: https://github.com/michaelschueler/scitools/releases/latest
github: https://github.com/michaelschueler
display: public
         protected
source: true
search: true

**SciTools** is a Fortran library of mathematical and numerical tools for
scientific applications. It builds a single library, `libscitools`; all
modules are prefixed `scitools_`. For installation and build options see the
[README](https://github.com/michaelschueler/scitools#readme).

This documentation describes the public interface of each module.

Modules
-------

| Area | Modules |
|------|---------|
| Basics | [[scitools_def]] (kinds, constants, units), [[scitools_utils]] (strings, I/O helpers, errors), [[scitools_debug]] (assertions) |
| Linear algebra | [[scitools_linalg]] (BLAS/LAPACK wrappers) |
| Iterative linear solvers | [[scitools_cg_solver]], [[scitools_bicgstab]], [[scitools_gmresm]], [[scitools_idrs]] |
| Quadrature | [[scitools_gausslegendre]], [[scitools_gauss_lobatto]], [[scitools_gregory]], [[scitools_filon]], [[scitools_quadpack]], [[scitools_quadrature]], [[scitools_lebedev_quad]], [[scitools_lebedev_weights]] |
| Interpolation and differentiation | [[scitools_bsplines]], [[scitools_vector_bsplines]], [[scitools_fornberg]] |
| ODEs and time evolution | [[scitools_rungekutta]], [[scitools_matrixexp]], [[scitools_evol]] |
| Special functions and polynomials | [[scitools_special]], [[scitools_legendre]], [[scitools_hermite]] |
| Root finding, sorting, random numbers | [[scitools_root]], [[scitools_sorting]], [[scitools_random]] |
| Timing | [[scitools_time]] |
| File I/O | [[scitools_binutils]], [[scitools_hdf5_utils]] (requires HDF5) |
| Parallelization | [[scitools_array1d_dist]] (requires MPI) |
| Physics helpers | [[scitools_laserpulse]], [[scitools_gausspulse]] |

Releases
--------

The current stable release can be
[downloaded on GitHub](https://github.com/michaelschueler/scitools/releases/latest).
All releases and the change log are listed on the
[releases page](https://github.com/michaelschueler/scitools/releases).
