# Third-party notices

SciTools is licensed under the Mozilla Public License 2.0 (see `LICENSE.md`).
Parts of SciTools are based on, or adapted from, the works listed below, which
remain under their original licenses. The full license texts are in the
`licenses/` directory.

| SciTools source | Original work | Authors | License | License text |
|-----------------|---------------|---------|---------|--------------|
| `src/basic/linalg.f90`, `src/basic/utils.F90`, `src/basic/sorting.f90`, `src/basic/debug.F90` (shape assertions), `src/math/special.F90` | [fortran-utils](https://github.com/certik/fortran-utils); utilities and sorting partly based on code by J. E. Pask (LLNL) | O. Čertík | MIT | `licenses/fortran-utils.txt` |
| `src/math/special.F90` (`SPHJ`, `SPHY`, `MSTA1`, `MSTA2`) | `specfun.f` (via SciPy and fortran-utils) | S. Zhang, J. Jin | Permission to use in programs | `licenses/specfun.txt` |
| `src/math/bsplines.f90` | [bspline-fortran](https://github.com/jacobwilliams/bspline-fortran), based on NIST CMLIB and SLATEC | J. Williams; NIST; C. de Boor | BSD-3-Clause; public domain | `licenses/bspline-fortran.txt`, `licenses/public-domain.txt` |
| `src/math/gausslegendre.F90` | [quadrature-fortran](https://github.com/jacobwilliams/quadrature-fortran), including SLATEC `DGAUS8` | J. Williams; SLATEC | BSD-3-Clause; public domain | `licenses/quadrature-fortran.txt`, `licenses/public-domain.txt` |
| `src/math/quadpack.f90` | QUADPACK | R. Piessens, E. de Doncker-Kapenger, C. Überhuber, D. Kahaner | Public domain | `licenses/public-domain.txt` |
| `src/math/hermite.f90`, `src/math/filon.F90` | Routines by J. Burkardt (`hermite` partly after S. Elhay and J. Kautsky) | J. Burkardt | MIT | `licenses/burkardt-mit.txt` |

<!-- TODO before release: confirm origin and license of the following -->
Still to be confirmed:

* `src/math/idrs.f90`: adapted from the IDR(s) code of M. B. van Gijzen
  (P. Sonneveld and M. B. van Gijzen, SIAM J. Sci. Comput. 31, 1035 (2008)).
* `src/math/gmresm.f90`: GMRES(m) with hookstep, apparently after A. P. Willis.
* Data tables in `src/math/lebedev/` (Lebedev-Laikov grids) and `src/math/lobatto/`.

Implemented by the SciTools authors from the literature (no third-party code):

* `src/math/matrixexp.F90`: scaling and squaring with Padé approximants,
  N. J. Higham, SIAM J. Matrix Anal. Appl. 26, 1179 (2005).
* `src/math/bicgstab.F90`: H. A. van der Vorst, SIAM J. Sci. Stat. Comput. 13, 631 (1992).
* `src/math/fornberg.f90`: B. Fornberg, Math. Comp. 51, 699 (1988).
