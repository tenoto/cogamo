#!/usr/bin/env python

import os 

datadir = '/Users/enoto/Dropbox/01_enoto/research/growth/logbook/220310_Cogamo_response/data/miwa/20220323/'

outdir = 'outdir/data_20220323'
cmd = 'rm -rf %s;mkdir -p %s' % (outdir,outdir)
os.system(cmd) 

## plot curve 
cmd = '$COGAMO_PATH/cogamo/cli/cgm_plot_curve.py '
cmd += '%s/032_20211215_15_137Cs.csv ' % datadir
cmd += '--tbin 1.0 '
cmd += '--tstart 2300.0 --tstop 3100.0'
print(cmd);#os.system(cmd)

cmd = '$COGAMO_PATH/cogamo/cli/cgm_plot_curve.py '
cmd += '%s/032_20211215_17_60Co.csv ' % datadir
print(cmd);os.system(cmd)

cmd = '$COGAMO_PATH/cogamo/cli/cgm_plot_curve.py '
cmd += '%s/032_20211215_16_Bgd.csv ' % datadir
print(cmd);os.system(cmd)

## plot spectrum
cmd = '$COGAMO_PATH/cogamo/cli/cgm_plot_adc_spectrum.py '
cmd += '%s/032_20211215_15_137Cs.csv ' % datadir
print(cmd);os.system(cmd)

cmd = '$COGAMO_PATH/cogamo/cli/cgm_plot_adc_spectrum.py '
cmd += '%s/032_20211215_17_60Co.csv ' % datadir
print(cmd);os.system(cmd)

cmd = '$COGAMO_PATH/cogamo/cli/cgm_plot_adc_spectrum.py '
cmd += '%s/032_20211215_16_Bgd.csv ' % datadir
print(cmd);os.system(cmd)

cmd = '$COGAMO_PATH/cogamo/cli/cgm_fit_phaspec_line.py '
cmd += '%s/032_20211215_15_137Cs.csv ' % datadir
print(cmd);os.system(cmd)

os.system('mv *.pdf *.yaml %s/' % outdir)