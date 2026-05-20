cmake_minimum_required(VERSION 3.10)
project("SCITOOLS" Fortran)
enable_language (Fortran)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)


set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# ~~ Commenly shared libs ~~
set(COMMON_LIBS "")

# Require MPI for this project:
if(mpi)
   message(STATUS "Building with MPI")
   find_package(MPI REQUIRED COMPONENTS Fortran)
   list(APPEND COMMON_LIBS MPI::MPI_Fortran)
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DMPI")
endif(mpi)

# ~~ Add HDF5 ~~
if (hdf5)
  message(STATUS "Building with HDF5")
  find_package(HDF5 REQUIRED COMPONENTS Fortran)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DWITHHDF5")
  include_directories(${HDF5_INCLUDE_DIRS})
  list(APPEND COMMON_LIBS ${HDF5_Fortran_LIBRARIES} ${HDF5_LIBRARIES})
endif (hdf5)

# ~~ Add BLAS/LAPACK ~~
if(lapackblas_libraries)
   list(APPEND COMMON_LIBS ${lapackblas_libraries})
   message(STATUS LAPACK/BLAS: ${lapackblas_libraries})
else(lapackblas_libraries)
   message("Warning: LAPACK/BLAS not specified!")
endif(lapackblas_libraries)

# ~~ Add modern quadrature ~~
if (modernquad)
  message(STATUS "Building with modern quadrature")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DMODERNQUAD")
endif (modernquad)

# ~~ CMake Configuration ~~
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ~~ Libraries ~~
link_directories(${CMAKE_LIBRARY_PATH})

# ~~ Paths and Subdirs ~~
include_directories(${CMAKE_INCLUDE_PATH})

add_subdirectory(src)


# ~~ Include Guard ~~
set(CONFIGURED_ONCE TRUE CACHE INTERNAL
    "A flag showing, that CMake has configured at least once.")
