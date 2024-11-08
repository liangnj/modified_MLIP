#!/bin/bash
#SBATCH -o job.%j.out
#SBATCH -e job.%j.err
#SBATCH -p C064M0256G
#SBATCH --qos=low
#SBATCH -J lammps
#SBATCH -N 4
#SBATCH --ntasks=256
#SBATCH --ntasks-per-node=64

module load mkl/latest
module load compiler/latest
module load mpi/latest

mpirun -np 256 /lustre/home/2001111732/software/lammps-MTP1/mlip-2-master/lammps-8Feb2023/src/lmp_intel_cpu_intelmpi -in graphene.in > graphene.out
