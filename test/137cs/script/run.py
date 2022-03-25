#!/usr/bin/env python

import os 

import cogamo.cogamo as cogamo

datadir = '/Users/enoto/Dropbox/01_enoto/research/growth/logbook/220310_Cogamo_response/data/miwa/20220323/'

file_bgd = '%s/032_20211215_16_Bgd.csv' % datadir
file_137cs = '%s/032_20211215_15_137Cs.csv' % datadir
file_60co = '%s/032_20211215_17_60Co.csv' % datadir

DICT_INITPAR_K40 = {'name':'K40','MeV':1.46083,'peak':142,'sigma':3,'area':2000,'c0':500,'c1':-3.0,'pha_min':100,'pha_max':164,'binning':2,'xlim':[100,164]}
DICT_INITPAR_TL208 = {'name':'Tl208','MeV':2.61453,'peak':250,'sigma':2,'area':520,'c0':70,'c1':-0.3,'pha_min':230,'pha_max':294,'binning':2,'xlim':[230,294],}
GAMMA_LINES = [DICT_INITPAR_K40,DICT_INITPAR_TL208]

evt_bgd = cogamo.EventData(file_bgd)
evt_bgd.set_outdir('out/bgd')
evt_bgd.extract_pha_spectrum(binning=2)
evt_bgd.prepare_energy_calibration(GAMMA_LINES)
evt_bgd.set_energy_series()
evt_bgd.set_time_series()
evt_bgd.plot_multi_curves()
evt_bgd.extract_energy_spectrum()
evt_bgd.extract_xspec_pha()

evt_137cs = cogamo.EventData(file_137cs)
evt_137cs.set_outdir('out/137cs')
evt_137cs.extract_pha_spectrum(binning=2)
evt_137cs.set_energy_series(pha2mev_c1=evt_bgd.param['pha2mev_c1'],pha2mev_c0=evt_bgd.param['pha2mev_c0'])
evt_137cs.set_time_series()
evt_137cs.plot_multi_curves()
evt_137cs.extract_energy_spectrum()
evt_137cs.extract_xspec_pha()

evt_60co = cogamo.EventData(file_60co)
evt_60co.set_outdir('out/60co')
evt_60co.extract_pha_spectrum(binning=2)
evt_60co.set_energy_series(pha2mev_c1=evt_bgd.param['pha2mev_c1'],pha2mev_c0=evt_bgd.param['pha2mev_c0'])
evt_60co.set_time_series()
evt_60co.plot_multi_curves()
evt_60co.extract_energy_spectrum()
evt_60co.extract_xspec_pha()
