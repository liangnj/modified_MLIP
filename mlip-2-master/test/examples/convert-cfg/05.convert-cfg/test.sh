#!/bin/bash

# Preamble, common for all examples
MLP_EXE=../../../../bin/mlp
TMP_DIR=./out
mkdir -p $TMP_DIR

$MLP_EXE convert-cfg input.bin.cfg out/POSCAR --input-format=bin --output-format=vasp-poscar > /dev/null
diff correct_POSCAR0 out/POSCAR0 1>&2
diff correct_POSCAR1 out/POSCAR1 1>&2
diff correct_POSCAR2 out/POSCAR2 1>&2
diff correct_POSCAR3 out/POSCAR3 1>&2

