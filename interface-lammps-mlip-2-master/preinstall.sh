#!/bin/sh

LAMMPS_PATH=$1

if [ ! -f "lib_mlip_interface.a" ]; then echo Error: copy lib_mlip_interface.a to this directory.; echo "       Create it by doing make libinterface in the MLIP root folder"; exit 1; fi
if [ ! -d "$LAMMPS_PATH/src" ]; then echo Error: run ./preinstall.sh \<path-to-lammps\>; exit 1; fi

curdir=$(pwd)
echo $curdir
mkdir -p $LAMMPS_PATH/lib/mlip
cp lib_mlip_interface.a $LAMMPS_PATH/lib/mlip
cp LAMMPS/Makefile.lammps.template $LAMMPS_PATH/lib/mlip/Makefile.lammps
cp -r LAMMPS/USER-MLIP/ $LAMMPS_PATH/src/
cp LAMMPS/Install.sh $LAMMPS_PATH/src/USER-MLIP/
cp LAMMPS/README $LAMMPS_PATH/src/USER-MLIP/
cd $LAMMPS_PATH/src
make no-user-mlip
make yes-user-mlip

