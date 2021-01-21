include ${GYROEP_ROOT}/shared/install/make.inc.${GYROEP_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

EXEC = locpargen

EXTRA_LIBS = ${GYROEP_ROOT}/shared/EXPRO/EXPRO_lib.a \
             ${GYROEP_ROOT}/shared/GEO/GEO_lib.a  \
             ${GYROEP_ROOT}/shared/math/math_lib.a 

OBJECTS = $(EXEC).o

.SUFFIXES : .o .f90

exec: $(OBJECTS)
	$(FC) $(FFLAGS) -o $(EXEC) $(OBJECTS) $(EXTRA_LIBS) $(LMATH)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -f *.o $(EXEC)
	rm -f ${GYROEP_ROOT}/modules/locpargen*.mod
