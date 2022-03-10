#!/bin/sh -f 

rm -rf out
mkdir out
cgm_plot_hkfile.py data/038_20211219.csv
cgm_plot_hkfile.py data/012_20210109.csv 
cgm_plot_hkfile.py data/014_20210109.csv

mv 038_20211219_hk.pdf out/
mv 012_20210109_hk.pdf out/
mv 014_20210109_hk.pdf out/