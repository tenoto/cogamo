#!/bin/sh -f

rm -rf  out 
mkdir ./out 

$COGAMO_PATH/cogamo/response/01_fill_eventlist_to_matrix.py

mv 01_matrix.pdf  01_matrix.txt 01_matrix_zoom.pdf out/

