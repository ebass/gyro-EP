include ${GYROEP_ROOT}/shared/install/make.inc.${GYROEP_PLATFORM}

export EXTRA_LIBS = \
	${GYROEP_ROOT}/shared/math/math_lib.a \
	${GYROEP_ROOT}/shared/EXPRO/EXPRO_lib.a \
	${GYROEP_ROOT}/shared/GEO/GEO_lib.a \
	${LMATH} \
	${NETCDF}

all:
	cd ${GYROEP_ROOT}/shared/math ; make
	cd ${GYROEP_ROOT}/shared/GEO ; make
	cd ${GYROEP_ROOT}/shared/EXPRO ; make
	cd locpargen ; make
	cd src ; make
	cd exprotools ; make

clean:
	cd locpargen ; make clean
	cd src ; make clean
	cd exprotools ; make clean
