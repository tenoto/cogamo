#!/usr/bin/env python

import cogamo.cogamo as cogamo

outdir = 'out'

evt = cogamo.EventData('../data/event/id014/014_20211219_14.csv')
#evt = cogamo.EventData('../data/event/id056/056_20211219_14.csv')
#evt = cogamo.EventData('../data/event/id056/056_20211219_00.csv') # low bgd, non detection 
outdir = '%s/product/id%s/%s/%s/%s/%s' % (outdir,evt.detid_str, evt.year, evt.month, evt.day, evt.hour_jst)

evt.set_outdir(outdir)

evt.extract_pha_spectrum(binning=2)
evt.prepare_energy_calibration()
evt.set_energy_series()
evt.set_time_series()
evt.plot_multi_curves()
evt.extract_energy_spectrum()
evt.extract_xspec_pha()

lc = evt.extract_curve()
evt.search_burst(lc,threshold_sigma=5.0)
if evt.numof_bst > 0:
	evt.analysis_bursts()

evt.show_summary()
evt.write_to_fitsfile()
evt.write_to_yamlfile()	
evt.pdfmerge()
print(evt.param)
