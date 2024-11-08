#/bin/bash

# Preamble, common for all examples
MLP_EXE=../../../bin/mlp
TMP_DIR=./out
mkdir -p $TMP_DIR

# Body:
$MLP_EXE calc-grade pot.mtp train.cfg test.cfg $TMP_DIR/graded.cfg --als-filename=$TMP_DIR/state.als

