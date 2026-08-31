cmake_minimum_required(VERSION 3.14)

project(sl3d
  LANGUAGES Fortran
  DESCRIPTION "Forward and adjoint modelling of glacial isostatic adjustment")

# Both sub-projects are independent programs sharing the same build settings.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# The original Makefiles always compiled with -fbounds-check; keep that as the
# default, but allow it to be switched off for production runs.
option(SL3D_BOUNDS_CHECK "Compile with runtime array bounds checking" ON)

include(GNUInstallDirs)

# Everything the build produces stays inside the build tree: executables and
# input data in build/run, libraries in build/lib, and .mod files alongside
# them.  Nothing is written back into the source directories.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

find_package(OpenMP REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_library(FFTW3_LIBRARY NAMES fftw3)

# netCDF-Fortran, used to read the ICE5G ice histories in data/ice5g.
# On Debian or Ubuntu this is the libnetcdff-dev package.
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
  pkg_check_modules(NETCDF_FORTRAN QUIET IMPORTED_TARGET netcdf-fortran)
endif()
if(NOT NETCDF_FORTRAN_FOUND)
  message(FATAL_ERROR
    "netCDF-Fortran not found. Install it (libnetcdff-dev on Debian/Ubuntu).")
endif()

# Common flags for the legacy sources: fixed-form .f files and some
# type-punned calls need gfortran's relaxed argument checking (>= 10).
add_library(sl3d_flags INTERFACE)
target_link_libraries(sl3d_flags INTERFACE
  OpenMP::OpenMP_Fortran
  LAPACK::LAPACK
  BLAS::BLAS)
if(FFTW3_LIBRARY)
  target_link_libraries(sl3d_flags INTERFACE ${FFTW3_LIBRARY})
endif()
target_link_libraries(sl3d_flags INTERFACE PkgConfig::NETCDF_FORTRAN)

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  target_compile_options(sl3d_flags INTERFACE
    -fbacktrace
    $<$<VERSION_GREATER_EQUAL:$<Fortran_COMPILER_VERSION>,10>:-fallow-argument-mismatch>
    $<$<BOOL:${SL3D_BOUNDS_CHECK}>:-fbounds-check>)
endif()

# Adds one program of <code> (gia1d or sl3d), built from src/<name>.f90 and
# linked against that code's core library.  Targets are named <code>_<name>
# because the two codes share several program names (merge_files, visc_test,
# file_test); for the same reason each code keeps its own directory inside
# build/run, so the programs themselves keep their familiar names.
function(sl3d_add_program code name)
  add_executable(${code}_${name} src/${name}.f90)
  target_link_libraries(${code}_${name} PRIVATE ${code}_core)
  set_target_properties(${code}_${name} PROPERTIES
    OUTPUT_NAME ${name}
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod/${name}
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/run/${code})
  install(TARGETS ${code}_${name}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/${code})
endfunction()

add_subdirectory(gia1d)
add_subdirectory(sl3d)

# The programs read their inputs from, and write their outputs to, the current
# working directory, so build/run is laid out as the place to run them from:
# the input data beside the executables, and results kept out of the source
# tree.
#
#   cd build/run && ./sl3d/sea_level_3d prem.simp.vm1
#
file(GLOB SL3D_INPUT_DATA CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/data/*)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/run)

# The ice histories are read only and rather larger, so link rather than copy.
foreach(ice ice5g ice6g ice7g)
  list(REMOVE_ITEM SL3D_INPUT_DATA ${CMAKE_CURRENT_SOURCE_DIR}/data/${ice})
  file(CREATE_LINK
    ${CMAKE_CURRENT_SOURCE_DIR}/data/${ice} ${CMAKE_BINARY_DIR}/run/${ice}
    SYMBOLIC COPY_ON_ERROR)
endforeach()

file(COPY ${SL3D_INPUT_DATA} DESTINATION ${CMAKE_BINARY_DIR}/run)
