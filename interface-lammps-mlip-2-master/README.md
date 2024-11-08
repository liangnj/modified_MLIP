# interface-lammps-mlip-v2

An interface between LAMMPS and MLIP.
MLIP, Machine Learning Interatomic Potentials, is a software distributed
by the group of A. Shapeev from Skoltech (Moscow).
MLIP can be obtained by following the instructions at
https://mlip.skoltech.ru/download/.
After you have obtained MLIP you can follow the instructions below to
make the MLIP-LAMMPS interaface.

Note: current version of the MLIP-LAMMPS interface works with 
LAMMPS-2022 or later. Please, download one of these versions of LAMMPS. 

## Installation

In order to install a MLIP plugin into LAMMPS:

0. Clone MLIP-LAMMPS interface

```git clone https://gitlab.com/ashapeev/interface-lammps-mlip-2.git .```

1. Clone the 2022 (or later) version of the LAMMPS code (we recommend using a previously uncompiled version of the code) 

```git clone -b stable https://github.com/lammps/lammps.git mylammps```
 
and place it in the same folder as `interface-lammps-mlip-2/`

2. Compile `lib_mlip_interface.a` in the MLIP folder

```make libinterface```

and place it in `interface-lammps-mlip-2/`

3. (optional) If additional LAMMPS packages are necessary they should be included at the end of the preinstall.sh file

```make yes-<package>```

e.g.

```make yes-manybody```

4. Compile LAMMPS  

```./install.sh <path-to-lammps> <lammps-target>```

e.g.

```./install.sh ../mylammps intel_cpu_intelmpi```

If everything goes fine, the lammps binary will be copied to this folder

### Combination of LAMMPS and MLIP settings that work:

* MLIP and LAMMPS with default compiler (serial version) `./install.sh ../mylammps serial` 
* MLIP and LAMMPS with default compiler (parallel version) `./install.sh ../mylammps mpi` 
* MLIP with the GNU compiler and LAMMPS with `./install.sh ../mylammps g++_serial` 
* MLIP with the GNU compiler and mpich and LAMMPS with `./install.sh ../mylammps g++_mpich` 
* MLIP with Intel compiler and intel MPI and LAMMPS with `./install.sh ../mylammps intel_cpu_intelmpi`
* If you'd like to build with Intel and no MPI then you need to find a way to link LAMMPS with MKL
