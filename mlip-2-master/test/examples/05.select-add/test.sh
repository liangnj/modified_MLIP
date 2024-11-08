#/bin/bash

# Preamble, common for all examples
MLP_EXE=../../../bin/mlp
TMP_DIR=./out
mkdir -p $TMP_DIR

# Body:
$MLP_EXE select-add pot.mtp train.cfg candidate.cfg $TMP_DIR/to-be-computed-on-dft.cfg --als-filename=$TMP_DIR/state.als --selected-filename=$TMP_DIR/selected.cfg
