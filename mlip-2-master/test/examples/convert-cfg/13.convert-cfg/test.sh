#!/bin/bash

# Preamble, common for all examples
MLP_EXE=../../../../bin/mlp
TMP_DIR=./out
mkdir -p $TMP_DIR

cp output.cfg out/output.cfg
$MLP_EXE convert-cfg OUTCAR out/output.cfg --input-format=vasp-outcar --append > /dev/null
diff correct_output.cfg out/output.cfg 1>&2
