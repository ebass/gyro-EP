include ${GYROEP_ROOT}/shared/install/make.inc.${GYROEP_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

LLIB=TRANSP_lib

OBJECTS = fTRANSP_GLOBALS.o \
          fTRANSP_INIT.o \
          fTRANSP_DO.o \
          fTRANSP_CLEANUP.o \
          rTRANSP_GLOBALS.o \
          rTRANSP_INIT.o \
          rTRANSP_DO.o \
          rTRANSP_CLEANUP.o 

.SUFFIXES : .o .f90

lib: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS) 

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:  
	rm -rf *.o $(LLIB).a
	rm -rf ${GYROEP_ROOT}/modules/ftransp*.mod
	rm -rf ${GYROEP_ROOT}/modules/rtransp*.mod
	                
