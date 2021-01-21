#--------------------------------------------
# Environment variable setup for gyro-EP
#--------------------------------------------
#!/bin/tcsh

setenv PATH ${PATH}:$GYROEP_ROOT/bin
setenv PATH ${PATH}:$GYROEP_ROOT/shared/bin

if ( $?PYTHONPATH ) then
 setenv PYTHONPATH ${PYTHONPATH}:$GYROEP_ROOT/shared/python
else
 setenv PYTHONPATH $GYROEP_ROOT/shared/python
endif

