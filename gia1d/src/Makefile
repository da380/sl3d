MOD=../mod/
LOBJ=../obj/
BIN=../bin/
MV=mv -f
F95=gfortran
DHOME=/home/da380



LIB= -L/usr/local/lapack-3.8.0 -llapack -lblas \
     -L$(DHOME)/raid/SHTOOLS/lib -lSHTOOLS2.8 \
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



all:    merge_files                                   \
        sea_level_1d                                  \
        sl_init_inv                                   \
        grace_visc_kern                               \
        ice_kern_1d                                   \
        visc_kern_1d


merge_files: $(BIN)merge_files $(OBJS)
$(BIN)merge_files: $(OBJS) $(LOBJ)merge_files.o
	$(F95) $(FFLAGS) $(LOBJ)merge_files.o $(OBJS) $(LIB) -o $(BIN)merge_files

$(LOBJ)merge_files.o: merge_files.f90 
	$(F95) $(FFLAGS) merge_files.f90 -c
	$(MV) merge_files.o $(LOBJ)

sea_level_1d: $(BIN)sea_level_1d $(OBJS)
$(BIN)sea_level_1d: $(OBJS) $(LOBJ)sea_level_1d.o
	$(F95) $(FFLAGS) $(LOBJ)sea_level_1d.o $(OBJS) $(LIB) -o $(BIN)sea_level_1d

$(LOBJ)sea_level_1d.o: sea_level_1d.f90 
	$(F95) $(FFLAGS) sea_level_1d.f90 -c
	$(MV) sea_level_1d.o $(LOBJ)

adj_test: $(BIN)adj_test $(OBJS)
$(BIN)adj_test: $(OBJS) $(LOBJ)adj_test.o
	$(F95) $(FFLAGS) $(LOBJ)adj_test.o $(OBJS) $(LIB) -o $(BIN)adj_test

$(LOBJ)adj_test.o: adj_test.f90 
	$(F95) $(FFLAGS) adj_test.f90 -c
	$(MV) adj_test.o $(LOBJ)

visc_test: $(BIN)visc_test $(OBJS)
$(BIN)visc_test: $(OBJS) $(LOBJ)visc_test.o
	$(F95) $(FFLAGS) $(LOBJ)visc_test.o $(OBJS) $(LIB) -o $(BIN)visc_test

$(LOBJ)visc_test.o: visc_test.f90 
	$(F95) $(FFLAGS) visc_test.f90 -c
	$(MV) visc_test.o $(LOBJ)

sl_init_inv: $(BIN)sl_init_inv $(OBJS)
$(BIN)sl_init_inv: $(OBJS) $(LOBJ)sl_init_inv.o
	$(F95) $(FFLAGS) $(LOBJ)sl_init_inv.o $(OBJS) $(LIB) -o $(BIN)sl_init_inv

$(LOBJ)sl_init_inv.o: sl_init_inv.f90 
	$(F95) $(FFLAGS) sl_init_inv.f90 -c
	$(MV) sl_init_inv.o $(LOBJ)

grace_visc_kern: $(BIN)grace_visc_kern $(OBJS)
$(BIN)grace_visc_kern: $(OBJS) $(LOBJ)grace_visc_kern.o
	$(F95) $(FFLAGS) $(LOBJ)grace_visc_kern.o $(OBJS) $(LIB) -o $(BIN)grace_visc_kern

$(LOBJ)grace_visc_kern.o: grace_visc_kern.f90 
	$(F95) $(FFLAGS) grace_visc_kern.f90 -c
	$(MV) grace_visc_kern.o $(LOBJ)

ice_kern_1d: $(BIN)ice_kern_1d $(OBJS)
$(BIN)ice_kern_1d: $(OBJS) $(LOBJ)ice_kern_1d.o
	$(F95) $(FFLAGS) $(LOBJ)ice_kern_1d.o $(OBJS) $(LIB) -o $(BIN)ice_kern_1d

$(LOBJ)ice_kern_1d.o: ice_kern_1d.f90 
	$(F95) $(FFLAGS) ice_kern_1d.f90 -c
	$(MV) ice_kern_1d.o $(LOBJ)

visc_kern_1d: $(BIN)visc_kern_1d $(OBJS)
$(BIN)visc_kern_1d: $(OBJS) $(LOBJ)visc_kern_1d.o
	$(F95) $(FFLAGS) $(LOBJ)visc_kern_1d.o $(OBJS) $(LIB) -o $(BIN)visc_kern_1d

$(LOBJ)visc_kern_1d.o: visc_kern_1d.f90 
	$(F95) $(FFLAGS) visc_kern_1d.f90 -c
	$(MV) visc_kern_1d.o $(LOBJ)



$(OBJS): $(SRCS)	
	$(F95)   $(FFLAGS)  -c $(SRCS)
	$(MV)  *.o $(LOBJ)
