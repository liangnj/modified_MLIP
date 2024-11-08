#/bin/bash

# Preamble, common for all examples
MLP_EXE=../../../bin/mlp
TMP_DIR=./out
mkdir -p $TMP_DIR

# Body:
$MLP_EXE relax mlip.ini --cfg-filename=to-relax.cfg --save-relaxed=$TMP_DIR/relaxed.cfg
