#!/bin/bash

# Preamble, common for all examples
MLP_EXE=../../../../bin/mlp
TMP_DIR=./out
mkdir -p $TMP_DIR

$MLP_EXE convert-cfg input.cfg out/POSCAR --output-format=vasp-poscar > /dev/null 
diff correctPOSCAR out/POSCAR 1>&2

