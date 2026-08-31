# gia1d -- GIA in 1D viscosity models (with 3D sensitivity kernels).
# Source list and ordering follow the original src/Makefile, which is now kept
# only in the git history.
# module_sphexp.f90 is deliberately excluded: it wraps SHTOOLS, which the
# original Makefile also left out of the object list, and nothing else uses it.

set(GIA1D_CORE_SRCS
  src/nrtype.f90
  src/nrutil.f90
  src/module_util.f90
  src/module_ice.f90
  src/module_config.f90
  src/module_fourier.f90
  src/module_model.f90
  src/gll_library.f90
  src/lagrange_poly.f90
  src/get_global.f90
  src/module_spline.f90
  src/rotmx2.f
  src/prott.f
  src/module_function.f90
  src/module_mat.f90
  src/module_gia.f90
  src/module_sl.f90
  src/module_kern.f90)

add_library(gia1d_core STATIC ${GIA1D_CORE_SRCS})
target_link_libraries(gia1d_core PUBLIC sl3d_flags)
set_target_properties(gia1d_core PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(gia1d_core PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/mod)

set(GIA1D_PROGRAMS
  merge_files
  sea_level_1d
  sl_init_inv
  grace_visc_kern
  ice_kern_1d
  visc_kern_1d
  visc_test
  file_test)

# Not built: src/adj_test.f90 and src/kernel_int.f90.  Neither was in the
# original Makefile's "all" target and both have fallen out of sync with the
# modules: kernel_int calls visc_3d, which only exists in sl3d, and adj_test
# calls forward_full_sl/stat_adj/adj_sl_ice with argument lists that no longer
# match module_sl.  Reviving them means updating the science, not the build.

foreach(prog IN LISTS GIA1D_PROGRAMS)
  sl3d_add_program(gia1d ${prog})
endforeach()
