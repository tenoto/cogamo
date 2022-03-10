#!/bin/bash 

echo "#################"
echo "#    COGAMO     #"
echo "#################"

export COGAMO_PATH=$(pwd)
export PYTHONPATH=$COGAMO_PATH:$PYTHONPATH
echo "COGAMO_PATH=" $COGAMO_PATH
export PATH=$COGAMO_PATH/cogamo/cli:$PATH

export COGAMO_SAMPLE_DATA_PATH="/Users/enoto/Dropbox/01_enoto/research/growth/data/20211130_otsuru_sample/cogamo_data_v20211130"
echo "COGAMO_SAMPLE_DATA_PATH=" $COGAMO_SAMPLE_DATA_PATH

export COGAMO_RESPONSE_DATA_PATH="/Users/enoto/Dropbox/01_enoto/research/growth/logbook/220310_Cogamo_response/data/gabriel/Reponse_files"
echo "COGAMO_RESPONSE_DATA_PATH=" $COGAMO_RESPONSE_DATA_PATH

