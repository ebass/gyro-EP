include ${GYROEP_ROOT}/shared/install/make.inc.${GYROEP_PLATFORM}

all:
	cd ${GYROEP_ROOT}/shared/math ; make
	cd ${GYROEP_ROOT}/shared/GEO ; make
	cd ${GYROEP_ROOT}/shared/EXPRO ; make
	cd ${GYROEP_ROOT}/shared/fluxfit ; make
	cd locpargen ; make
	cd src ; make
	cd exprotools ; make

clean:
	cd locpargen ; make clean
	cd src ; make clean
	cd exprotools ; make clean
