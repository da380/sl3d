MOD=../mod/
LOBJ=../obj/
BIN=../bin/
MV=mv -f
F95=gfortran
DHOME=/home/da380



LIB= -llapack -L$(DHOME)/raid/SHTOOLS/lib -lSHTOOLS2.8 \
     -lfftw3 -lm -s -m64


FFLAGS =  -I../mod -J../mod -O5          \
          -I$(DHOME)/raid/SHTOOLS/modules \
          -fbounds-check                  \
	  -fbacktrace                     \
          -fopenmp

SRCS =  nrtype.f90                                    \
	nrutil.f90                                    \
	module_util.f90                               \
	module_fourier.f90                            \
        module_model.f90                              \
	gll_library.f90                               \
	lagrange_poly.f90                             \
        get_global.f90                                \
        module_spline.f90                             \
        rotmx2.f                                      \
        prott.f	                                      \
	rsple.f                                       \
        rspln.f                                       \
        module_function.f90                           \
        module_mat.f90                                \
        module_gia.f90                                \
        module_sl.f90                                 \
        module_kern.f90



OBJS =  $(LOBJ)nrtype.o                               \
	$(LOBJ)nrutil.o                               \
	$(LOBJ)module_util.o                          \
        $(LOBJ)module_fourier.o                       \
        $(LOBJ)module_model.o                         \
	$(LOBJ)gll_library.o                          \
	$(LOBJ)lagrange_poly.o                        \
        $(LOBJ)get_global.o                           \
        $(LOBJ)module_spline.o                        \
        $(LOBJ)rotmx2.o                               \
        $(LOBJ)prott.o                                \
        $(LOBJ)rsple.o                                \
        $(LOBJ)rspln.o                                \
        $(LOBJ)module_function.o                      \
        $(LOBJ)module_mat.o                           \
        $(LOBJ)module_gia.o                           \
        $(LOBJ)module_sl.o                            \
        $(LOBJ)module_kern.o



MODS =  $(MOD)nrtype.mod                              \
	$(MOD)nrutil.mod                              \
	$(MOD)module_util.mod                         \
        $(MOD)module_fourier.mod                      \
        $(MOD)module_model.mod                        \
        $(MOD)module_function.mod                     \
        $(MOD)module_mat.mod                          \
        $(MOD)module_spline.mod                       \
        $(MOD)module_gia.mod                          \
        $(MOD)module_sl.mod                           \
        $(MOD)module_kern.mod



all:    visco_test                                    \
        merge_files                                   \
        sea_level_3d                                  \
        kernel_int                                    \
        read_kl                                       \
        read_kl_2                                     \
        sl_visc_inv                                   \
        grace_visc_inv                                \
        ice_kern_3d                                   \
        visc_kern_3d                                  \
        rsl_inversion                                 \
        read_rob                                      \
        visc_corr


visco_test: $(BIN)visco_test $(OBJS)
$(BIN)visco_test: $(OBJS) $(LOBJ)visco_test.o
	$(F95) $(FFLAGS) $(LOBJ)visco_test.o $(OBJS) $(LIB) -o $(BIN)visco_test

$(LOBJ)visco_test.o: visco_test.f90 
	$(F95) $(FFLAGS) visco_test.f90 -c
	$(MV) visco_test.o $(LOBJ)

merge_files: $(BIN)merge_files $(OBJS)
$(BIN)merge_files: $(OBJS) $(LOBJ)merge_files.o
	$(F95) $(FFLAGS) $(LOBJ)merge_files.o $(OBJS) $(LIB) -o $(BIN)merge_files

$(LOBJ)merge_files.o: merge_files.f90 
	$(F95) $(FFLAGS) merge_files.f90 -c
	$(MV) merge_files.o $(LOBJ)

sea_level_3d: $(BIN)sea_level_3d $(OBJS)
$(BIN)sea_level_3d: $(OBJS) $(LOBJ)sea_level_3d.o
	$(F95) $(FFLAGS) $(LOBJ)sea_level_3d.o $(OBJS) $(LIB) -o $(BIN)sea_level_3d

$(LOBJ)sea_level_3d.o: sea_level_3d.f90 
	$(F95) $(FFLAGS) sea_level_3d.f90 -c
	$(MV) sea_level_3d.o $(LOBJ)

kernel_int: $(BIN)kernel_int $(OBJS)
$(BIN)kernel_int: $(OBJS) $(LOBJ)kernel_int.o
	$(F95) $(FFLAGS) $(LOBJ)kernel_int.o $(OBJS) $(LIB) -o $(BIN)kernel_int

$(LOBJ)kernel_int.o: kernel_int.f90 
	$(F95) $(FFLAGS) kernel_int.f90 -c
	$(MV) kernel_int.o $(LOBJ)

read_kl: $(BIN)read_kl $(OBJS)
$(BIN)read_kl: $(OBJS) $(LOBJ)read_kl.o
	$(F95) $(FFLAGS) $(LOBJ)read_kl.o $(OBJS) $(LIB) -o $(BIN)read_kl

$(LOBJ)read_kl.o: read_kl.f90 
	$(F95) $(FFLAGS) read_kl.f90 -c
	$(MV) read_kl.o $(LOBJ)

read_kl_2: $(BIN)read_kl_2 $(OBJS)
$(BIN)read_kl_2: $(OBJS) $(LOBJ)read_kl_2.o
	$(F95) $(FFLAGS) $(LOBJ)read_kl_2.o $(OBJS) $(LIB) -o $(BIN)read_kl_2

$(LOBJ)read_kl_2.o: read_kl_2.f90 
	$(F95) $(FFLAGS) read_kl_2.f90 -c
	$(MV) read_kl_2.o $(LOBJ)

sl_visc_inv: $(BIN)sl_visc_inv $(OBJS)
$(BIN)sl_visc_inv: $(OBJS) $(LOBJ)sl_visc_inv.o
	$(F95) $(FFLAGS) $(LOBJ)sl_visc_inv.o $(OBJS) $(LIB) -o $(BIN)sl_visc_inv

$(LOBJ)sl_visc_inv.o: sl_visc_inv.f90 
	$(F95) $(FFLAGS) sl_visc_inv.f90 -c
	$(MV) sl_visc_inv.o $(LOBJ)

grace_visc_inv: $(BIN)grace_visc_inv $(OBJS)
$(BIN)grace_visc_inv: $(OBJS) $(LOBJ)grace_visc_inv.o
	$(F95) $(FFLAGS) $(LOBJ)grace_visc_inv.o $(OBJS) $(LIB) -o $(BIN)grace_visc_inv

$(LOBJ)grace_visc_inv.o: grace_visc_inv.f90 
	$(F95) $(FFLAGS) grace_visc_inv.f90 -c
	$(MV) grace_visc_inv.o $(LOBJ)

ice_kern_3d: $(BIN)ice_kern_3d $(OBJS)
$(BIN)ice_kern_3d: $(OBJS) $(LOBJ)ice_kern_3d.o
	$(F95) $(FFLAGS) $(LOBJ)ice_kern_3d.o $(OBJS) $(LIB) -o $(BIN)ice_kern_3d

$(LOBJ)ice_kern_3d.o: ice_kern_3d.f90 
	$(F95) $(FFLAGS) ice_kern_3d.f90 -c
	$(MV) ice_kern_3d.o $(LOBJ)

visc_kern_3d: $(BIN)visc_kern_3d $(OBJS)
$(BIN)visc_kern_3d: $(OBJS) $(LOBJ)visc_kern_3d.o
	$(F95) $(FFLAGS) $(LOBJ)visc_kern_3d.o $(OBJS) $(LIB) -o $(BIN)visc_kern_3d

$(LOBJ)visc_kern_3d.o: visc_kern_3d.f90 
	$(F95) $(FFLAGS) visc_kern_3d.f90 -c
	$(MV) visc_kern_3d.o $(LOBJ)

rsl_inversion: $(BIN)rsl_inversion $(OBJS)
$(BIN)rsl_inversion: $(OBJS) $(LOBJ)rsl_inversion.o
	$(F95) $(FFLAGS) $(LOBJ)rsl_inversion.o $(OBJS) $(LIB) -o $(BIN)rsl_inversion

$(LOBJ)rsl_inversion.o: rsl_inversion.f90 
	$(F95) $(FFLAGS) rsl_inversion.f90 -c
	$(MV) rsl_inversion.o $(LOBJ)

read_rob: $(BIN)read_rob $(OBJS)
$(BIN)read_rob: $(OBJS) $(LOBJ)read_rob.o
	$(F95) $(FFLAGS) $(LOBJ)read_rob.o $(OBJS) $(LIB) -o $(BIN)read_rob

$(LOBJ)read_rob.o: read_rob.f90 
	$(F95) $(FFLAGS) read_rob.f90 -c
	$(MV) read_rob.o $(LOBJ)

visc_corr: $(BIN)visc_corr $(OBJS)
$(BIN)visc_corr: $(OBJS) $(LOBJ)visc_corr.o
	$(F95) $(FFLAGS) $(LOBJ)visc_corr.o $(OBJS) $(LIB) -o $(BIN)visc_corr

$(LOBJ)visc_corr.o: visc_corr.f90 
	$(F95) $(FFLAGS) visc_corr.f90 -c
	$(MV) visc_corr.o $(LOBJ)

$(OBJS): $(SRCS)	
	$(F95)   $(FFLAGS)  -c $(SRCS)
	$(MV)  *.o $(LOBJ)
