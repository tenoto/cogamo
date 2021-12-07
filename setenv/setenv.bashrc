#!/bin/bash 

echo "#################"
echo "#    COGAMO     #"
echo "#################"

export COGAMO_PATH=$(pwd)
export PYTHONPATH=$COGAMO_PATH:$PYTHONPATH

echo "COGAMO_PATH:" $COGAMO_PATH