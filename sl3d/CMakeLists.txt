# sl3d -- GIA in 3D viscosity models via a pseudospectral method.
# Source list and ordering follow the original src/Makefile, which is now kept
# only in the git history.
# module_sphexp.f90 (SHTOOLS wrappers) and module_kern_stress.f90 are excluded,
# as they were in the original object list.

set(SL3D_CORE_SRCS
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
  src/rsple.f90
  src/rspln.f90
  src/module_function.f90
  src/module_mat.f90
  src/module_gia.f90
  src/module_sl.f90
  src/module_kern.f90)

add_library(sl3d_core STATIC ${SL3D_CORE_SRCS})
target_link_libraries(sl3d_core PUBLIC sl3d_flags)
set_target_properties(sl3d_core PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)
target_include_directories(sl3d_core PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/mod)

set(SL3D_PROGRAMS
  visco_test
  merge_files
  sea_level_3d
  kernel_int
  read_kl
  read_kl_2
  sl_visc_inv
  grace_visc_inv
  ice_kern_3d
  visc_kern_3d
  rsl_inversion
  read_rob
  visc_corr
  adj_test
  visc_test
  extract_sl
  file_test)

foreach(prog IN LISTS SL3D_PROGRAMS)
  sl3d_add_program(sl3d ${prog})
endforeach()
