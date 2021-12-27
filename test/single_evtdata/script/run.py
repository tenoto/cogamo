#!/usr/bin/env python

import numpy as np 
import cogamo.cogamo as cogamo

outdir = 'out'

evt = cogamo.EventData('/Users/enoto/Dropbox/01_enoto/research/growth/data/20211130_otsuru_sample/cogamo_data_v20211130/014_20210109_06.csv')
#evt = cogamo.EventData('../data/012_20210108_17.csv')
outdir = '%s/product/id%s/%s/%s/%s/%s' % (outdir,evt.detid_str, evt.year, evt.month, evt.day, evt.hour_jst)
evt.set_outdir(outdir)
evt.extract_pha_spectrum(binning=2)
evt.prepare_energy_calibration()
evt.set_energy_series()
evt.set_time_series()
evt.extract_energy_spectrum()
lc = evt.extract_curve()
evt.search_burst(lc,threshold_sigma=4.0)
if evt.numof_bst > 0:
	evt.analysis_bursts()
evt.show_summary()
evt.write_to_fitsfile()
evt.write_to_yamlfile()	
evt.pdfmerge()
print(evt.param)
