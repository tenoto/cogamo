#!/bin/sh 

cgm_plot_phasepc.py data/20220411/032_20220411_17_Bgd.csv --tstart 2720 --tstop 2730 --output_pdf 032_20220411_17_bgd_spec_2720to2730.pdf --flag_rate 1 --ymin 0.01 --ymax 50 
cgm_plot_phasepc.py data/20220411/032_20220411_17_Bgd.csv --tstart 2260 --tstop 2340 --output_pdf 032_20220411_17_bgd_spec_2260to2340.pdf --flag_rate 1 --ymin 0.01 --ymax 50 

cgm_plot_curve.py data/20220411/032_20220411_17_Bgd.csv --phamin 0 --tstart 2100 --tstop 2800 --tbin 1
cgm_plot_curve.py data/20220411/032_20220411_17_Bgd.csv --phamin 0 --tstart 2100 --tstop 2800 --tbin 1 --phamin 100 --output_pdf 032_20220411_17_bgd_curve_phamin100.pdf 
