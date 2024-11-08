#!/bin/bash

# Preamble, common for all examples
MLP_EXE=../../../../bin/mlp
TMP_DIR=./out
mkdir -p $TMP_DIR

$MLP_EXE convert-cfg OUTCAR out/POSCAR --input-format=vasp-outcar --output-format=vasp-poscar > /dev/null
diff correct_POSCAR out/POSCAR 1>&2

