# COMPILATION AND INSTALLATION PROCEDURE

## QUICK OVERVIEW
MLIP has a serial and parallel version. Building the serial version of MLIP
requires a modern C++ compiler (supporting c++11), parallel version requires
an MPI C++ compiler with c++11, as well as C and FORTRAN compilers.

MLIP uses the BLAS library. MLIP can work without this library using only
the embedded version of BLAS (although you are advised to use the
high-performance libraries, e.g., Intel MKL or OpenBLAS).

MLIP build system consists of two step: configuration and compilation.
The configuration is performed by running:
```bash
./configure [OPTIONS]
```

Detailed descriptions of configuration options can be viewed by:
```bash
./configure --help
```

After configuration, the executables can be built with make command:
```bash
make <target-name>
``` 
The target will be created in bin/ or lib/

For information about the build executables, execute
```bash
make help
```

## DETAILED INSTRUCTIONS

1. Configure the MLIP package build
First, start with executing the configuration script
```bash
./configure
```
The configuration attempts to set the correct values for various
system-dependent variables for compilation. The values related to BLAS are
autodetected in the case if Intel MKL is installed on the system.
If not, then the "embedded" BLAS (supplied with MLIP) will be used.
You can force to use the embedded BLAS with the `--blas=embedded` option:
```bash
./configure --blas=embedded
```
By default, the parallel version of the MLIP is compiled. Compilation of serial
version can be done by specifying the `--no-mpi` option to configure':
```bash
./configure --no-mpi
```

2. Build the MLIP executables
You can compile and build the MLIP executable by command:
```bash
make mlp
```

3. Install OpenBLAS (optional)
You may want to download and install the OpenBLAS library:
```bash
git clone https://github.com/xianyi/OpenBLAS.git
make -C OpenBLAS 
make PREFIX=./ -C OpenBLAS install
```
(You, of course, can download it directly from 
https://github.com/xianyi/OpenBLAS).
You should then configure MLIP with
```bash
./configure --blas=openblas --blas-root=<OpenBLAS path>
```

## BUILDING MLIP LIBRARY FOR INTERFACES WITH OTHER CODES (E.G. LAMMPS)

```bash
make libinterface
```
This creates `lib/lib_mlip_interface.a` which should be used when building
other packages embedding MLIP.
The interface with MLIP is avaiable
[here](https://gitlab.com/ashapeev/interface-lammps-mlip-2/)
under a GNU Public License.

## TESTING THE MLIP BUILD
The build can be tested by running:
```bash 
make test
```
This launches a number of tests, including test of LAMMPS if it is present in `bin/`
(i.e., if `bin/lmp_serial` and `bin/lmp_mpi` exist).

