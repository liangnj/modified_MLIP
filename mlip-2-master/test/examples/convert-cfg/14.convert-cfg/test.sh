#!/bin/bash

# Preamble, common for all examples
MLP_EXE=../../../../bin/mlp
TMP_DIR=./out
mkdir -p $TMP_DIR

$MLP_EXE convert-cfg OUTCAR_relax out/output.cfg --input-format=vasp-outcar --output-format=txt > /dev/null 2> /dev/null
diff correct_output.cfg out/output.cfg 1>&2
