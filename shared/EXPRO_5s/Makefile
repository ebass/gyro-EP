include ${GYROEP_ROOT}/shared/install/make.inc.${GYROEP_PLATFORM}

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

LLIB=EXPRO_lib

export EXPRO_EXTRA_LIBS = \
       ${GYROEP_ROOT}/shared/math/math_lib.a

OBJECTS = EXPRO_interface.o \
          EXPRO_globals.o \
          EXPRO_alloc.o \
          EXPRO_palloc.o \
          EXPRO_alloc_control.o \
          EXPRO_read.o \
          EXPRO_pread.o \
          EXPRO_read_driver.o \
          EXPRO_compute_derived.o \
          EXPRO_write_derived.o \
          EXPRO_write_original.o

.SUFFIXES : .o .f90 

all: lib

lib: $(OBJECTS)
	$(ARCH) $(LLIB).a ${EXPRO_EXTRA_LIBS} $(OBJECTS)

.f90.o :
	$(FC) $(FMATH) $(FFLAGS) -c $<

clean:
	rm -f *.o $(LLIB).a ${GYROEP_ROOT}/modules/expro_*.mod

