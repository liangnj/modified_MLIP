#/bin/bash

# Preamble, common for all examples
MLP_EXE=../../../bin/mlp
TMP_DIR=./out
mkdir -p $TMP_DIR

# Body:
$MLP_EXE train 06.mtp train.cfg --trained-pot-name=$TMP_DIR/pot.mtp --max-iter=100

