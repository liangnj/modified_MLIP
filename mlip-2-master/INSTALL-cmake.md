# COMPILATION AND INSTALLATION PROCEDURE

## QUICK OVERVIEW
MLIP has a serial and parallel version. Building the serial version of MLIP
requires a modern C++ compiler (supporting c++11), parallel version requires
an MPI C++ compiler with c++11, as well as C and FORTRAN compilers.

MLIP uses the BLAS library. MLIP can work without this library using only
the embedded version of BLAS (although you are advised to use the
high-performance libraries, e.g., Intel MKL or OpenBLAS).

## Build MLIP with CMake
The following text assumes some familiarity with CMake and focuses on using
the command line tool cmake and what settings are supported for building
MLIP.

**NOTE:**
MLIP currently requires that CMake version 3.10 or later is available.

### Getting started
MLIP build system is a two-step process:

- First, in a new directory using either the command-line utility (cmake)
or the text-mode UI utility ccmake or the graphical utility cmake-gui,
generate a build environment. E.g., one can use the command-line version of
the CMake with no customization as,
   ```bash
   cd mlip-2                # change to the MLIP distribution directory
   mkdir build              # create a build directory
   cd build
   cmake ..                 # configuration
   ```
   This process generates build-files for the default build command. During the
   configuration step, CMake will try to detect whether support for MPI and BLAS
   are available and enable the corresponding configuration settings.

- Second step is compiling the code, and linking of objects, libraries, and
executables.
   ```bash
   make                     # compilation
   ```
   The make command launches the compilation, and, if successful, will ultimately
   produce a library `lib_mlip_interface.a` and the MLIP executable `mlp` inside
   the build folder.

   You can speed this process by a parallel compilation with `make -j N` (with `N`
   being the number of concurrently executed tasks when your machine has multiple
   CPU cores). For example,
   ```bash
   make -j4                 # compilation
   ```

   After compilation, you may optionally install the MLIP executable and created
   libraries into your system with:
   ```bash
   make install             # optional, copy files into installation location
   ```

## DETAILED INSTRUCTIONS

You can specify some options that affect the build process in the MLIP
executable/interface.

**NOTE:**

A required argument is a directory, which contains the `CMakeLists.txt` file.
For example,
```bash
cmake /path/to/mlip-2 [options]   # build from any dir
cmake .. [options]                # build from mlip-2/build
```

The CMake command-line options useful for building MLIP are,

```bash
-DCMAKE_INSTALL_PREFIX=path  # where to install MLIP executable, lib if desired

-DCMAKE_BUILD_TYPE=type      # type = Debug, Release, RelWithDebInfo, MinSizeRel, default is 'RelWithDebInfo' if nothing is specified

-DWITH_MPI=value             # ON or OFF, default is ON if CMake finds MPI, else OFF
-DWITH_SELFTEST=value        # ON or OFF, default is ON (Enable 'self-testing' implementation)
-DWITH_LOSSLESS_CFG=value    # ON or OFF, default is OFF (Enable 'lossless configuration' file writing)
-DWITH_SYSTEM_BLAS=value     # ON or OFF, default is ON (Enable using system built BLAS library (with a C-style interface))
-DBLAS_ROOT=path             # where to look for system/user built BLAS library (with a C-style interface)
```

**NOTE: CMAKE_BUILD_TYPE OPTION**

`CMAKE_BUILD_TYPE` specifies what build type (configuration) will be
built in this build tree. Selecting the build type causes CMake to
set many per-config properties and variables. For example, in a build
tree configured to build `Release,`

```bash
cmake ../cmake -DCMAKE_BUILD_TYPE=Release
```

CMake will see having RELEASE settings get added to the compilers FLAGS.
Thus when configuring MLIP in this mode, optimization flags are used.

Or in a build tree configured to build type Debug,

```bash
cmake ../cmake -DCMAKE_BUILD_TYPE=Debug
```

CMake will see to having DEBUG settings get added to the compilers FLAGS.

If no building type is specified, it defaults to a `RelWithDebInfo` (short
for Release With Debug Information) build.

**NOTE: WITH_MPI OPTION**

The default option `WITH_MPI=ON`, is to find MPI and build a parallel version
of MLIP. By setting this option to `OFF`, one can avoid the search for MPI
and force to build a serial version.

**NOTE: WITH_SYSTEM_BLAS OPTION**

The default option `WITH_SYSTEM_BLAS=ON`, is to use the system installed `BLAS`
library and its C interface. The priority when using the `Intel compilers` is
the `MKL` package if it exists, where its environment variables are set via
(`MKLROOT`).

Next the build process searches for `OpenBLAS` installation, and next, it looks
for other libraries like `Atlas` or `LAPACK` and checks if it can use them. In
case no library is found or the found library and the header files are not usable
by this process, or the required header files are missing, it would build the
embedded library and use it.

In the case of `WITH_SYSTEM_BLAS=OFF` value, there is no search for the system
built library and it builds its own internal library from the embedded source
files and use it.

To use the specific user/system installed library, you can set a `BLAS_ROOT`
variable to a directory that contains a `BLAS` installation. Setting this
option precedes any other system search. Usually the `BLAS_ROOT` path is the
root to the installation directory and should contain the `lib` and `include`
sub-directories.

### Choice of compiler and compile/link options

An important option to gain performance is the choice of compilers (and
compiler flags). Some vendor compilers like
[Intel compilers](https://software.intel.com/content/www/us/en/develop/tools/compilers.html)
can produce faster code than open-source compilers like GNU.

The CMake command-line options useful for compiler and compile/link options
are,
```bash
-DCMAKE_CXX_COMPILER=name            # A name of C++ compiler, e.g., icpc, g++
-DCMAKE_C_COMPILER=name              # A name of C compiler, e.g., icc, gcc
-DCMAKE_Fortran_COMPILER=name        # A name of Fortran compiler, e.g., ifort, gfortran

-DCMAKE_CXX_FLAGS=string             # C++ compiler flags
-DCMAKE_C_FLAGS=string               # C compiler flags
-DCMAKE_Fortran_FLAGS=string         # Fortran compiler flags
```

For example, for building with Intel Compilers, within build from mlip-2/build
directoy
```bash
cmake ../cmake -DCMAKE_CXX_COMPILER=icpc -DCMAKE_C_COMPILER=icc -DCMAKE_Fortran_COMPILER=ifort
```

or, for building with the GNU Compilers (latest version), within build from
mlip-2/build directoy
```bash
cmake ../cmake -DCMAKE_CXX_COMPILER=g++-10 -DCMAKE_C_COMPILER=gcc-10 -DCMAKE_Fortran_COMPILER=gfortran-10
```

or, for building with the specific MPI wrapper compilers, within build from mlip-2/build
directoy
```bash
cmake ../cmake -DCMAKE_CXX_COMPILER=/path/to/mpi/bin/mpic++ -DCMAKE_C_COMPILER=/path/to/mpi/bin/mpicc -DCMAKE_Fortran_COMPILER=/path/to/mpi/bin/mpifort
```

**NOTE:**

After the build process you can compile only the MLIP executable by,
```bash
make mlp
```

You can also compile and build only the MLIP interface to use with other
codes like LAMMPS by,
```bash
make libinterface
```

## BUILDING MLIP LIBRARY FOR INTERFACES WITH OTHER CODES (E.G. LAMMPS)

Compiling the code, and linking of objects, libraries, and executables by,
```bash
make
```
if successful, will ultimately produce a library `lib_mlip_interface.a` and
the MLIP executable `mlp` inside the build folder.

The `lib_mlip_interface.a` which should be used when building other packages
embedding MLIP.

The interface with MLIP is avaiable
[here](https://gitlab.com/ashapeev/interface-lammps-mlip-2/)
under a GNU Public License.


## Install OpenBLAS (optional)
You may want to download and install the OpenBLAS library:
```bash
git clone https://github.com/xianyi/OpenBLAS.git
make -C OpenBLAS
make PREFIX=./ -C OpenBLAS install
```
(You, of course, can download it directly from
[OpenBLAS](https://github.com/xianyi/OpenBLAS).)


You can then build MLIP with
```bash
   mkdir build                          # create a build directory
   cd build
   cmake .. -DBLAS_ROOT=<OpenBLAS path> # configuration
   make -j4
```

## TESTING THE MLIP BUILD
The build can be tested by running:
```bash
make test
```
This launches a number of tests.
