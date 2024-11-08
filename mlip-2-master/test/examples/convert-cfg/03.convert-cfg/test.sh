#!/bin/bash

# Preamble, common for all examples
MLP_EXE=../../../../bin/mlp
TMP_DIR=./out
mkdir -p $TMP_DIR

$MLP_EXE convert-cfg input.cfg out/output.bin.cfg --output-format=bin > /dev/null 
diff correct_output.bin.cfg out/output.bin.cfg 1>&2

