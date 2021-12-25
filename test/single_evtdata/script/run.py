#!/usr/bin/env python

import cogamo.cogamo as cogamo

outdir = 'out'

##################################################
# Default value 
##################################################

PHA_SPECTRUM_NBINS = 2**10

ENERGY_SPECTRUM_NBINS = 2**10
ENERGY_SPECTRUM_MIN = 0.0
ENERGY_SPECTRUM_MAX = 12.0

PEAK_MEV_K40 = 1.46083
PEAK_MEV_TL208 = 2.61453

DICT_INITPAR_TL208 = {'peak':236,'sigma':7,'area':2651,'c0':798.0,'c1':-3,'MeV':PEAK_MEV_TL208,'pha_min':200,'pha_max':284,'nbins':83,'xlim':[200,284],'name':'Tl208'}
DICT_INITPAR_K40 = {'peak':132,'sigma':5,'area':18025,'c0':3731,'c1':-21.0,'MeV':PEAK_MEV_K40,'pha_min':100,'pha_max':164,'nbins':64,'xlim':[100,164],'name':'K40'}

DICT_EXTRACT_CURVE = {'tbin':5.0,'tstart':0.0,'tstop':3600.0,'energy_min':3.0,'energy_max':None,'xlim':[0.0,3600.0]}
DICT_SEARCH_BURST = DICT_EXTRACT_CURVE
DICT_SEARCH_BURST['burst_sigma'] = 4.0
DICT_FIT_BURST_CURVE = DICT_EXTRACT_CURVE
DICT_FIT_BURST_CURVE['fit_nsigma'] = 8

##################################################
# Default value 
##################################################

evt = cogamo.EventData('../data/012_20210108_17.csv')
outdir = '%s/product/id%s/%s/%s/%s/%s' % (outdir,evt.detid_str, evt.year, evt.month, evt.day, evt.hour_jst)
evt.set_outdir(outdir)

evt.extract_pha_spectrum()
dict_par_Tl208 = evt.fit_pha_spectrum_line(DICT_INITPAR_TL208)
dict_par_K40 = evt.fit_pha_spectrum_line(DICT_INITPAR_K40)
evt.set_energy_calibration_curve(dict_par_K40,dict_par_Tl208)

evt.set_energy_series(
	pha2mev_c0=evt.pha2mev_c0,
	pha2mev_c1=evt.pha2mev_c1)
evt.set_time_series()

evt.extract_energy_spectrum()
evt.extract_curve()

evt.search_burst()

exit()

# if burst was detected
par = evt.fit_burst_curve()		

evt.get_burst_duration(par,tbin=1.0,tstart=par['fit_xmin'],tstop=par['fit_xmax'],
	energy_min=3.0,energy_max=None,linear_tbin_normalization=5.0)

evt.show_summary()

evt.write_to_fitsfile()
evt.write_to_yamlfile()




