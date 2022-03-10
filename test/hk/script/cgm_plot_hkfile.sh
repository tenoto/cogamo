#!/bin/sh -f 

rm -rf out
mkdir out
cgm_plot_hkfile.py ../data/log/id056/056_20211219.csv

mv 056_20211219_hk.pdf out/

