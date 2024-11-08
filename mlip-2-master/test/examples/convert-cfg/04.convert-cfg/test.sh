#!/bin/bash

# Preamble, common for all examples
MLP_EXE=../../../../bin/mlp
TMP_DIR=./out
mkdir -p $TMP_DIR

$MLP_EXE convert-cfg input.bin.cfg out/output.cfg --input-format=bin > /dev/null
diff correct_output.cfg out/output.cfg 1>&2

