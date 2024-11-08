#!/bin/bash

# Preamble, common for all examples
MLP_EXE=../../../../bin/mlp
TMP_DIR=./out
mkdir -p $TMP_DIR

$MLP_EXE convert-cfg sample_input.cfg out/lammps.inp --output-format=lammps-datafile --input-format=txt > /dev/null 2> /dev/null
diff correct_lammps.inp out/lammps.inp 1>&2

