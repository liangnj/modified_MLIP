#!/bin/bash

# Preamble, common for all examples
MLP_EXE=../../../../bin/mlp
TMP_DIR=./out
mkdir -p $TMP_DIR

$MLP_EXE convert-cfg input.cfg out/output.cfg --last > /dev/null
diff correct_output.cfg out/output.cfg 1>&2

