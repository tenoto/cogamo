# -*- coding: utf-8 -*-

import os 
import re
import sys
import yaml
import glob 
import numpy as np 
import pandas as pd 

from astropy.time import Time
from astropy.io import fits

from datetime import datetime, timedelta, timezone
tz_tokyo = timezone(timedelta(hours=+9), 'Asia/Tokyo')
tz_utc = timezone(timedelta(hours=0), 'UTC')

import matplotlib.pylab as plt 
import matplotlib.gridspec as gridspec

from scipy.stats import norm
from iminuit import Minuit
from probfit import Chi2Regression

##################################################
# Default value 
##################################################

PHA_SPECTRUM_NBINS = 2**10

ENERGY_SPECTRUM_NBINS = 2**10
ENERGY_SPECTRUM_MIN = 0.0
ENERGY_SPECTRUM_MAX = 12.0

PEAK_MEV_K40 = 1.46083
PEAK_MEV_TL208 = 2.61453

DICT_INITPAR_TL208 = {'peak':236,'sigma':7,'area':2651,'c0':798.0,'c1':-3,'MeV':PEAK_MEV_TL208,'pha_min':190,'pha_max':284,'nbins':64,'xlim':[190,284],'name':'Tl208'}
DICT_INITPAR_K40 = {'peak':132,'sigma':5,'area':18025,'c0':3731,'c1':-21.0,'MeV':PEAK_MEV_K40,'pha_min':100,'pha_max':164,'nbins':64,'xlim':[100,164],'name':'K40'}

DICT_EXTRACT_CURVE = {'tbin':8.0,'tstart':0.0,'tstop':3600.0,'energy_min':3.0,'energy_max':None,'xlim':[0.0,3600.0]}
DICT_SEARCH_BURST = DICT_EXTRACT_CURVE
DICT_SEARCH_BURST['threshold_sigma'] = 4.0
DICT_FIT_BURST_CURVE = DICT_EXTRACT_CURVE
DICT_FIT_BURST_CURVE['fit_nsigma'] = 8

##################################################
# General function and classes 
##################################################

def get_norm_probability(threshold,mean=0.0,standard_deviation=1.0):
	return 1.0-norm.cdf(threshold,mean,standard_deviation)

def plot_xydata(x,y,outpdf,yerr=None,model_x=None,model_y=None,
		xlabel='X title',ylabel='Y title',title='Title',
		flag_xlog=False,flag_ylog=False,xlim=None):

	fig, ax = plt.subplots(1,1, figsize=(11.69,8.27)) # A4 size, inich unit 
	fontsize = 18 	

	if yerr is not None:
		plt.errorbar(x,y,yerr=yerr,marker='o',ls='')
	else:
		plt.plot(x,y,marker='o',ls='')

	if (model_x is not None) and (model_y is not None):
		plt.plot(model_x,model_y,marker='',ls='--')

	plt.xlabel(xlabel, fontsize=fontsize)
	plt.ylabel(ylabel, fontsize=fontsize)
	plt.title(title,fontsize=fontsize)
	if flag_xlog: plt.xscale('log')				
	if flag_ylog: plt.yscale('log')
	if xlim!=None: plt.xlim(xlim)
	
	ax.minorticks_on()
	ax.grid(True)
	ax.grid(axis='both',which='major', linestyle='--', color='#000000')
	ax.grid(axis='both',which='minor', linestyle='--')	
	ax.tick_params(axis="both", which='major', direction='in', length=5)
	ax.tick_params(axis="both", which='minor', direction='in', length=3)

	plt.tick_params(labelsize=fontsize)
	plt.rcParams["font.family"] = "serif"
	plt.rcParams["mathtext.fontset"] = "dejavuserif"	
	plt.tight_layout()

	plt.savefig(outpdf)	

def plot_histogram(hist_x,hist_y,outpdf,hist_yerr=None,
		xlabel='X title',ylabel='Y title',title='Title',
		flag_xlog=False,flag_ylog=False,xlim=None,
		burst_mask=None):
	fig, ax = plt.subplots(1,1, figsize=(11.69,8.27)) # A4 size, inich unit 
	fontsize = 18 	

	if hist_yerr is not None:
		plt.errorbar(hist_x,hist_y,marker='',drawstyle='steps-mid',color='k')		
		plt.errorbar(hist_x,hist_y,yerr=hist_yerr,marker='',drawstyle='steps-mid',color='k')
	else:
		plt.errorbar(hist_x,hist_y,marker='',drawstyle='steps-mid',color='k')

	if burst_mask is not None:
		if hist_yerr is not None:
			plt.errorbar(hist_x[burst_mask],hist_y[burst_mask],marker='',drawstyle='steps-mid',color='r')		
			plt.errorbar(hist_x[burst_mask],hist_y[burst_mask],yerr=hist_yerr[burst_mask],marker='',drawstyle='steps-mid',color='r')
		else:
			plt.errorbar(hist_x[burst_mask],hist_y[burst_mask],marker='',drawstyle='steps-mid',color='r')		

	plt.xlabel(xlabel, fontsize=fontsize)
	plt.ylabel(ylabel, fontsize=fontsize)
	plt.title(title,fontsize=fontsize)
	if flag_xlog: plt.xscale('log')				
	if flag_ylog: plt.yscale('log')
	if xlim!=None: plt.xlim(xlim)
	
	ax.minorticks_on()
	ax.grid(True)
	ax.grid(axis='both',which='major', linestyle='--', color='#000000')
	ax.grid(axis='both',which='minor', linestyle='--')	
	ax.tick_params(axis="both", which='major', direction='in', length=5)
	ax.tick_params(axis="both", which='minor', direction='in', length=3)

	plt.tick_params(labelsize=fontsize)
	plt.rcParams["font.family"] = "serif"
	plt.rcParams["mathtext.fontset"] = "dejavuserif"	
	plt.tight_layout()

	plt.savefig(outpdf)	

def plot_fit_residual(
		hist_x,hist_y,model_x,model_y,outpdf,
		hist_xerr=None,hist_yerr=None,
		xlabel='X title',ylabel='Y title',title='Title',
		flag_xlog=False,flag_ylog=False,flag_hist=False,
		xlim=None,
		axvline_values=None,axvline_legends=None,
		legend_text='',legend_loc='lower left'):

	fig, axs = plt.subplots(2,1,figsize=(11.69,8.27), # A4 size, inich unit 
		sharex=True,gridspec_kw={'hspace':0},tight_layout=True)
	gs = gridspec.GridSpec(2,1,height_ratios=[3,1])
	gs.update(hspace=0.0)

	fontsize = 18
	fontsize_legend = 16
	plt.tight_layout()
	plt.tick_params(labelsize=fontsize)
	plt.rcParams["font.family"] = "serif"
	plt.rcParams["mathtext.fontset"] = "dejavuserif"	

	mask = hist_yerr > 0.0

	axs[0] = plt.subplot(gs[0,0])
	if flag_hist:
		axs[0].errorbar(hist_x,hist_y,
			xerr=hist_xerr,yerr=hist_yerr,
			marker=None,color='k',drawstyle='steps-mid')
	else:
		axs[0].errorbar(hist_x[mask],hist_y[mask],
			xerr=hist_xerr[mask],yerr=hist_yerr[mask],
			marker='o',linestyle='None',color='k')	
	axs[0].plot(model_x,model_y,c='red',drawstyle='steps-mid',linewidth=2)
	axs[0].set_ylabel(ylabel, fontsize=fontsize)
	axs[0].set_xlim(xlim)
	axs[0].set_title(title, fontsize=fontsize)	
	axs[0].get_xaxis().set_visible(False)

	axs[1] = plt.subplot(gs[1])
	if hist_yerr is not None:
		if flag_hist:
			axs[1].errorbar(
				hist_x[mask],(hist_y-model_y)[mask]/hist_yerr[mask],
				xerr=hist_xerr[mask],yerr=1.0,
				marker=None,linestyle='None',color='k')	
		else:
			axs[1].errorbar(
				hist_x[mask],(hist_y-model_y)[mask]/hist_yerr[mask],
				yerr=1.0,
				marker='o',linestyle='None',color='k')			
		axs[1].axhline(y=0,ls='--',color='r')
		axs[1].set_ylabel('(data-model)/error', fontsize=fontsize)
	axs[1].set_xlabel(xlabel, fontsize=fontsize)
	axs[1].set_xlim(xlim)

	for ax in axs:
		ax.minorticks_on()
		#ax.grid(True)
		#ax.grid(axis='both',which='major', linestyle='--', color='#000000')
		#ax.grid(axis='both',which='minor', linestyle='--')	
		ax.tick_params(axis="both", which='major', direction='in', length=5,labelsize=fontsize)
		ax.tick_params(axis="both", which='minor', direction='in', length=3,labelsize=fontsize)
		for value in axvline_values:
			if axvline_legends is None:
				ax.axvline(value,ls='--')				
			else:
				index = axvline_values.index(value)
				ax.axvline(value,ls='--',label=axvline_legends[index])
			
	legend = axs[0].legend(title=legend_text,loc=legend_loc,fontsize=fontsize_legend)
	legend.get_title().set_fontsize(fontsize_legend)
	#axs[1].legend(title="(data-model)/error",loc="lower right")

	fig.align_ylabels(axs)
	plt.savefig(outpdf)	

def model_gauss(x, peak, sigma, area):
    return area * np.exp(-0.5*(x-peak)**2/sigma**2)/(np.sqrt(2*np.pi)*sigma) 

def model_gauss_linear(x, peak, sigma, area, c0=0.0, c1=0.0):
    return area * np.exp(-0.5*(x-peak)**2/sigma**2)/(np.sqrt(2*np.pi)*sigma) + c0 + c1 * x

def model_linear(x, c0=0.0, c1=0.0):
    return c0 + c1 * x

def fit_gauss(x,y,par_init,error=None,fit_nsigma=3):

	if error is not None:
		chi2reg = Chi2Regression(model_gauss,x,y,error=error)
	else:
		chi2reg = Chi2Regression(model_gauss,x,y)	

	fit = Minuit(chi2reg, 
		peak=par_init['peak'],sigma=par_init['sigma'],area=par_init['area'],
		limit_peak=(0,None),limit_sigma=(0,None),limit_area=(0,None))
	fit.migrad()
	fit.minos() 
	#fit.print_param()

	peak = fit.values[0]
	sigma = fit.values[1]
	fit_xmin = peak - fit_nsigma * sigma 
	fit_xmax = peak + fit_nsigma * sigma 

	flag = np.logical_and(x >= fit_xmin, x <= fit_xmax)
	if error is not None:
		chi2reg = Chi2Regression(model_gauss,x[flag],y[flag],error=error[flag])
	else:
		chi2reg = Chi2Regression(model_gauss,x[flag],y[flag])

	fit = Minuit(chi2reg, 
		peak=par_init['peak'],sigma=par_init['sigma'],area=par_init['area'])
	fit.migrad()
	fit.minos() 
	fit.print_param()

	par = {
		'peak':fit.values[0],
		'sigma':fit.values[1],
		'area':fit.values[2],
		'c0':0.0,
		'c1':0.0,
		'peak_err':fit.errors[0],
		'sigma_err':fit.errors[1],
		'area_err':fit.errors[2],
		'c0_err':0.0,
		'c1_err':0.0,
		'fit_xmin':fit_xmin,
		'fit_xmax':fit_xmax,
		'fit_nsigma':fit_nsigma}
	return par

def fit_gauss_linear(x,y,par_init,error=None,fit_nsigma=3):

	if error is not None:
		chi2reg = Chi2Regression(model_gauss_linear,x,y,error=error)
	else:
		chi2reg = Chi2Regression(model_gauss_linear,x,y)	

	fit = Minuit(chi2reg, 
		peak=par_init['peak'],sigma=par_init['sigma'],area=par_init['area'], 
		c0=par_init['c0'], c1=par_init['c1'],
		limit_peak=(0,None),limit_sigma=(0,None),limit_area=(0,None))
	fit.migrad()
	fit.minos() 
	#fit.print_param()

	peak = fit.values[0]
	sigma = fit.values[1]
	fit_xmin = peak - fit_nsigma * sigma 
	fit_xmax = peak + fit_nsigma * sigma 

	flag = np.logical_and(x >= fit_xmin, x <= fit_xmax)
	if error is not None:
		chi2reg = Chi2Regression(model_gauss_linear,x[flag],y[flag],error=error[flag])
	else:
		chi2reg = Chi2Regression(model_gauss_linear,x[flag],y[flag])

	fit = Minuit(chi2reg, 
		peak=par_init['peak'],sigma=par_init['sigma'],area=par_init['area'], 
		c0=par_init['c0'], c1=par_init['c1'])
	fit.migrad()
	fit.minos() 
	fit.print_param()

	par = {
		'peak':fit.values[0],
		'sigma':fit.values[1],
		'area':fit.values[2],
		'c0':fit.values[3],
		'c1':fit.values[4],
		'peak_err':fit.errors[0],
		'sigma_err':fit.errors[1],
		'area_err':fit.errors[2],
		'c0_err':fit.errors[3],
		'c1_err':fit.errors[4],
		'fit_xmin':fit_xmin,
		'fit_xmax':fit_xmax,
		'fit_nsigma':fit_nsigma}
	return par

def get_energy_resolution(peak,sigma):
	# assuming the offset is zero
	return 2.35*sigma/peak 

class Hist1D(object):
	"""Represents 1D histogram. 
	"""
	def __init__(self, nbins, xlow, xhigh):
		self.nbins = nbins
		self.xlow  = xlow
		self.xhigh = xhigh
		self.y, edges = np.histogram([], bins=nbins, range=(xlow, xhigh))
		self.x = (edges[:-1] + edges[1:]) / 2.

		self.xerr = (-edges[:-1] + edges[1:]) / 2.

	def fill(self, arr):
		y, edges = np.histogram(arr, bins=self.nbins, range=(self.xlow, self.xhigh))
		self.y += y
		self.yerr = np.sqrt(self.y)

	@property
	def data(self):
		return self.x, self.y

##################################################
# Cogamo 
##################################################

class Archive(object):
	def __init__(self,parameter_yamlfile):
		self.param = yaml.load(open(parameter_yamlfile),
			Loader=yaml.FullLoader)		
		print("[Archive] %s is generatd." % self.param['archive_name'])

		self.dict = {
			'detid':[],
			'time':[],
			'process': [],
			'bstcand': [],
			'bstlclink': [],
			'bstlcfile': [],	
			'bstdistlink': [],
			'bstdistfile': [],	
			'bstlc_mean':[],
			'bstlc_sigma':[],
			'bstlc_area':[],	
			'bstalert_link': [],
			'bstalert_file': [],
			'lclink': [],
			'lcfile': [],
			'phalink': [],
			'phafile': [],
			'csvfile':[],
			'csvpath':[],			
			'csvlink':[]
			}

	def set_csvfiles(self):
		sys.stdout.write('[archive] {} \n'.format(sys._getframe().f_code.co_name))

		### Find CSV file and make archives. 
		csv_filelst = sorted(glob.glob('%s/**/*.csv' % self.param['datadir'],
			recursive=True))
		for file_path in csv_filelst:
			self.add(file_path)
		print(csv_filelst)	

	def add(self,file_path):
		print("[Archive %s] add %s" % (self.param['archive_name'],file_path))

		filename = os.path.basename(file_path)	
		basename, ext = os.path.splitext(filename)
		if ext != '.csv':
			print("Warning: skip the file %s" % csvfile)
			return -1
		detid, yyyymmdd, hour = basename.split("_")
		#print(detid,yyyymmdd,hour)

		year = yyyymmdd[0:4]
		month = yyyymmdd[4:6]
		day = yyyymmdd[6:8]		
		str_time = '%04d-%02d-%02dT%02d:00:00' % (int(year),int(month),int(day),int(hour))

		self.dict['detid'].append(detid)
		self.dict['time'].append(str_time)
		self.dict['process'].append('--')
		self.dict['bstcand'].append('--')
		self.dict['bstlclink'].append('--')
		self.dict['bstlcfile'].append('--')
		self.dict['bstdistlink'].append('--')
		self.dict['bstdistfile'].append('--')
		self.dict['bstlc_mean'].append('--')
		self.dict['bstlc_sigma'].append('--')
		self.dict['bstlc_area'].append('--')				
		self.dict['csvlink'].append('<a href="%s">%s</a>' % (file_path,filename))
		self.dict['csvpath'].append(file_path)		
		self.dict['csvfile'].append(filename)
		self.dict['bstalert_link'].append('--')
		self.dict['bstalert_file'].append('--')
		self.dict['lclink'].append('--')
		self.dict['lcfile'].append('--')			
		self.dict['phalink'].append('--')
		self.dict['phafile'].append('--')		
		return 0 

	def convert_to_dataframe(self):
		self.df = pd.DataFrame.from_dict(self.dict, orient='index').T
		#self.df = self.df.shift()[1:] # shift index starting from 0 to 1	

	def write(self):
		if not os.path.exists(self.param['outdir']):
			cmd = 'mkdir -p %s' % self.param['outdir']
			print(cmd);os.system(cmd)

		cmd = 'rm -f %s/%s.{csv,html}' % (self.param['outdir'],self.param['archive_name'])
		print(cmd);os.system(cmd)
	
		self.df.to_csv('%s/%s.csv' % (self.param['outdir'],self.param['archive_name']))

		self.df.drop(['csvpath','csvfile','lcfile','phafile','bstlcfile','bstdistfile','bstalert_file'],axis=1).to_html('%s/%s.html' % (self.param['outdir'],self.param['archive_name']), render_links=True, escape=False)

	def process(self,index):
		print("[Archive %s] process index of %s" % (self.param['archive_name'],index))

		csvpath = self.df.iloc[index]['csvpath']

		evt = EventData(csvpath)
		if evt.filetype != 'rawcsv':
			sys.stdout.write('Input file is not the rawcsv file format.')
			return 0		

		outdir = '%s/product/id%s/%s/%s/%s/%s' % (self.param['outdir'],evt.detid_str, evt.year, evt.month, evt.day, evt.hour_jst)
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

		# if burst was detected
		par = evt.fit_burst_curve()		

		evt.get_burst_duration(par,tbin=1.0,tstart=par['fit_xmin'],tstop=par['fit_xmax'],
			energy_min=3.0,energy_max=None,linear_tbin_normalization=5.0)

		evt.write_to_fitsfile()
		evt.write_to_yamlfile()

class EventData():
	"""Represents an EventData of a Cogamo detector.
	(raw csv, processed csv, or fits file formats)

	Raw CSV file (original file recorded by a cogamo detector):
		-- The event data file recorded all the individual radiation events.  
		-- The file name is [det_id]_[yyyymmdd]_[hour].csv (e.g., 012_20210108_17.csv)
		-- The time is defined in JST.
	The event data file include the following columns:
		1. minute
		2. sec
		3. 1/10000 sec
		4. ADC channel (pha: pulse height amplitude) [0-1023]
	Example: 
		0,0,77,60
		0,0,122,23
		0,0,166,41
	"""	
	def __init__(self, filepath):
		self.filetype = None		
		self.filepath = filepath
		self.filename = os.path.basename(self.filepath)
		self.basename = os.path.splitext(self.filename)[0]

		self.param = {}

		self.set_filetype()
		self.open_file()

	def set_filetype(self):
		if re.fullmatch(r'\d{3}_\d{8}_\d{2}.csv', self.filename):
			self.filetype = 'rawcsv'
			self.detid_str, self.yyyymmdd_jst, self.hour_jst = self.basename.split("_")		
			self.year = self.yyyymmdd_jst[0:4]
			self.month = self.yyyymmdd_jst[4:6]
			self.day = self.yyyymmdd_jst[6:8]			
		else:
			sys.stdout.write("[error] filetype error...")
			return -1

	def open_file(self):
		if not os.path.exists(self.filepath):
			raise FileNotFoundError("{} not found".format(self.filepath))
		try:
			if self.filetype == 'rawcsv':
				self.df = pd.read_csv(self.filepath, index_col=False, 
					names=['minute','sec','decisec','pha'],
					dtype={'minute':np.uintc,'sec':np.uintc,'decisec':np.uint16,'pha':np.uint16})
				self.nevents = len(self.df)
			else:
				sys.stdout.write("[error] filetype error...")
				return -1
		except OSError as e:
			raise

	def set_outdir(self,outdir,flag_overwrite=True):
		self.outdir = outdir 
		sys.stdout.write('----- {} -----\n'.format(sys._getframe().f_code.co_name))

		#outdir = 'out/%s' % evt.basename
		if flag_overwrite:
			cmd = 'rm -rf %s; mkdir -p %s' % (self.outdir,self.outdir)
			print(cmd);os.system(cmd)
		else:
			if os.path.exists(outdir):
				sys.stdout.write('outdir already existed.')
				return -1
			else:
				cmd = 'mkdir -p %s' % (self.outdir)
				print(cmd);os.system(cmd)				

	def set_time_series(self):		
		"""
		the standard unix time does not have enough accuracy below 1 second.
		Sub-second time stamp is handled by another column
		"""
		sys.stdout.write('----- {} -----\n'.format(sys._getframe().f_code.co_name))
		
		str_time = '%04d-%02d-%02dT%02d:' % (int(self.year),int(self.month),int(self.day),int(self.hour_jst))
		time_series_str  = np.char.array(np.full(self.nevents, str_time)) + np.char.mod('%02d:',self.df['minute']) + np.char.mod('%02d',self.df['sec']) + np.char.mod('.%04d',self.df['decisec']) 
		time_series_jst = Time(time_series_str, format='isot', scale='utc', precision=5) 
		self.time_series_utc = time_series_jst - timedelta(hours=+9)
		self.df['unixtime'] = self.time_series_utc.to_value('unix',subfmt='decimal')
		self.df['unixtime'] = self.df['unixtime'].astype(np.float64)

		time_offset_str = str_time + '00:00'
		time_offset_jst = Time(time_offset_str, format='isot', scale='utc', precision=5) 
		self.time_offset_utc = time_offset_jst - timedelta(hours=+9)		
		self.unixtime_offset = self.time_offset_utc.to_value('unix',subfmt='decimal')

	def set_energy_series(self,pha2mev_c0,pha2mev_c1):		
		""" with rand  
		"""
		sys.stdout.write('----- {} -----\n'.format(sys._getframe().f_code.co_name))

		rand_series = np.random.random_sample(len(self.df['pha']))
		rand_pha = self.df['pha'].astype(np.float32) + rand_series - 0.5 
		self.df['energy_mev'] = pha2mev_c1 * rand_pha + pha2mev_c0

	def extract_pha_spectrum(self,
			nbins=PHA_SPECTRUM_NBINS,pha_min=0.0,pha_max=PHA_SPECTRUM_NBINS,
			xlim=[1,PHA_SPECTRUM_NBINS]):
		"""Extract spectrum of pha
		nbins: pha bin size
		"""
		sys.stdout.write('----- {} -----\n'.format(sys._getframe().f_code.co_name))

		phaspec = PhaSpectrum(self.df['pha'],nbins=nbins,pha_min=pha_min,pha_max=pha_max)
		outpdf = '%s/%s_phaspec.pdf' % (self.outdir,self.basename)
		phaspec.plot(outpdf=outpdf,title=self.basename,xlim=xlim)
		return phaspec

	def extract_energy_spectrum(self,
			nbins=ENERGY_SPECTRUM_NBINS,
			energy_min=ENERGY_SPECTRUM_MIN,
			energy_max=ENERGY_SPECTRUM_MAX,
			xlim=[ENERGY_SPECTRUM_MIN,ENERGY_SPECTRUM_MAX]):
		"""Extract spectrum of pha
		nbins: pha bin size
		"""
		sys.stdout.write('----- {} -----\n'.format(sys._getframe().f_code.co_name))

		enespec = EnergySpectrum(self.df['energy_mev'],
			nbins=nbins,energy_min=energy_min,energy_max=energy_max)
		outpdf = '%s/%s_energyspec.pdf' % (self.outdir,self.basename)
		enespec.plot(outpdf=outpdf,title=self.basename,xlim=xlim)
		return enespec

	def fit_pha_spectrum_line(self,dict_initpar):
		"""Extract spectrum of pha
		nbins: pha bin size
		"""
		sys.stdout.write('----- {} -----\n'.format(sys._getframe().f_code.co_name))

		phaspec = PhaSpectrum(self.df['pha'],nbins=dict_initpar['nbins'],
			pha_min=dict_initpar['pha_min'],pha_max=dict_initpar['pha_max'])
		outpdf = '%s/%s_%s.pdf' % (self.outdir,self.basename,dict_initpar['name'])

		par = phaspec.fit_gauss_linear(dict_initpar,flag_error=True)
		par['name'] = dict_initpar['name']
		par['MeV'] = dict_initpar['MeV']		

		model_x = phaspec.hist.x
		model_y = np.array([model_gauss_linear(x,peak=par['peak'],sigma=par['sigma'],area=par['area'],c0=par['c0'],c1=par['c1']) for x in model_x])	

		legend_text = '%s at %.1f MeV (fit range: %d sigma)\n' % (par['name'],par['MeV'],par['fit_nsigma'])
		legend_text += 'peak=%.1f+/-%.1f ch\n' % (par['peak'],par['peak_err'])
		legend_text += 'sigma=%.1f+/-%.1f ch\n' % (par['sigma'],par['sigma_err'])		
		legend_text += 'area=%.1f+/-%.1f counts\n' % (par['area'],par['area_err'])
		legend_text += 'c0=%.1f+/-%.1f counts\n' % (par['c0'],par['c0_err'])		
		legend_text += 'c1=%.1f+/-%.1f counts\n' % (par['c1'],par['c1_err'])			
		legend_text += 'resolution=%.1f %%' % (100.0*get_energy_resolution(par['peak'],par['sigma']))

		plot_fit_residual(
			phaspec.hist.x,phaspec.hist.y,
			model_x,model_y,			
			outpdf=outpdf,
			hist_xerr=phaspec.hist.xerr,hist_yerr=phaspec.hist.yerr,
			xlabel='Channel',ylabel='Counts/bin',title=self.basename,
			xlim=[par['peak']-3.5*par['sigma'],par['peak']+3.5*par['sigma']],
			axvline_values=[par['fit_xmin'],par['fit_xmax']],
			legend_text=legend_text)
		return par

	def set_energy_calibration_curve(self,dict_par_K40,dict_par_Tl208):
		sys.stdout.write('----- {} -----\n'.format(sys._getframe().f_code.co_name))

		for par in [dict_par_K40,dict_par_Tl208]:
			for strtype in ['','_err']:
				self.param['%s_peak%s' % (par['name'],strtype)] = par['peak%s' % strtype]	
				self.param['%s_sigma%s' % (par['name'],strtype)] = par['sigma%s' % strtype]
				self.param['%s_area%s' % (par['name'],strtype)] = par['area%s' % strtype]
				self.param['%s_c0%s' % (par['name'],strtype)] = par['c0%s' % strtype]	
				self.param['%s_c1%s' % (par['name'],strtype)] = par['c1%s' % strtype]						
			self.param['%s_name' % par['name']] = par['name']
			self.param['%s_MeV' % par['name']] = par['MeV']								

		outpdf = '%s/%s_energycal.pdf' % (self.outdir,self.basename)
		self.pha2mev_c0, self.pha2mev_c1 = self.get_pha2mev_param(
			mev_array=np.array([PEAK_MEV_K40,PEAK_MEV_TL208]),
			pha_array=np.array([self.param['K40_peak'],self.param['Tl208_peak']]),
			pha_err_array=np.array([self.param['K40_peak_err'],self.param['Tl208_peak_err']]),
			title=self.basename,outpdf=outpdf)
		self.param['pha2mev_c0'] = self.pha2mev_c0
		self.param['pha2mev_c1'] = self.pha2mev_c1

	def get_pha2mev_param(self,mev_array,pha_array,pha_err_array,title,outpdf):
		sys.stdout.write('----- {} -----\n'.format(sys._getframe().f_code.co_name))

		mev_diff = np.diff(mev_array)[0]
		pha_diff = np.diff(pha_array)[0]
		slope_init = pha_diff / mev_diff
		chi2reg = Chi2Regression(model_linear,mev_array,pha_array,error=pha_err_array)
		fit = Minuit(chi2reg,c0=0.0,c1=slope_init)
		fit.migrad()
		fit.minos() 
		fit.print_param()	
		mev2pha_c0 = fit.values[0]	
		mev2pha_c1 = fit.values[1]	

		model_x = np.linspace(0,11,1600)
		model_y = np.array(model_linear(model_x,mev2pha_c0,mev2pha_c1))
		plot_xydata(mev_array,pha_array,yerr=pha_err_array,
			model_x=model_x, model_y=model_y,
			outpdf=outpdf,
			xlabel='Energy (MeV)',ylabel='PHA (channel)',title=title)

		pha2mev_c0 = - mev2pha_c0 / mev2pha_c1
		pha2mev_c1 = 1 / mev2pha_c1
		return pha2mev_c0, pha2mev_c1

	def get_energy_mask(self,energy_min,energy_max):

		if energy_min == None and energy_max == None:
			message = "no pha selection."		
			suffix = 'mev_all' 
			mask = np.full(len(self.df), True)
		elif energy_min != None and energy_max == None:
			message = "%.1f MeV <= Energy" % energy_min
			suffix = 'mev_%s_xx' % (str(energy_min).replace('.','p'))								
			mask = (self.df['energy_mev'] >= energy_min)
		elif energy_min == None and energy_max != None:
			message = "Energy <= %.1f MeV" % energy_max
			suffix = 'mev_xx_%s' % (str(energy_max).replace('.','p'))													
			mask = (self.df['energy_mev'] <= energy_max)
		elif energy_min != None and energy_max != None:
			message = "%.1f MeV <= Energy <= %.1f MeV" % (energy_min,energy_max)
			suffix = 'mev_%s_%s' % (str(energy_min).replace('.','p'),str(energy_max).replace('.','p'))			
			mask = np.logical_and((self.df['energy_mev'] >= energy_min),(self.df['energy_mev'] <= energy_max))
		sys.stdout.write("%s\n" % message)

		return mask, message, suffix 		

	def extract_curve(self,
			tbin=DICT_EXTRACT_CURVE['tbin'],
			tstart=DICT_EXTRACT_CURVE['tstart'],
			tstop=DICT_EXTRACT_CURVE['tstop'],
			energy_min=DICT_EXTRACT_CURVE['energy_min'],
			energy_max=DICT_EXTRACT_CURVE['energy_max'],
			xlim=DICT_EXTRACT_CURVE['xlim']):
		""" Energy in MeV unit
		"""
		sys.stdout.write('----- {} -----\n'.format(sys._getframe().f_code.co_name))

		mask, message, suffix = self.get_energy_mask(
			energy_min=energy_min,energy_max=energy_max)

		lc = LightCurve(
			np.array(self.df[mask]['unixtime']),
			float(self.unixtime_offset),
			tbin=tbin,tstart=tstart,tstop=tstop)
		title = '%s (%s)' % (self.basename, message)
		outpdf = '%s/%s_lc_%s.pdf' % (self.outdir,self.basename,suffix)
		lc.plot(outpdf=outpdf,title=title,xlim=xlim)		

	def search_burst(self,
			threshold_sigma=DICT_SEARCH_BURST['threshold_sigma'],
			tbin=DICT_SEARCH_BURST['tbin'],
			tstart=DICT_SEARCH_BURST['tstart'],
			tstop=DICT_SEARCH_BURST['tstop'],
			energy_min=DICT_SEARCH_BURST['energy_min'],
			energy_max=DICT_SEARCH_BURST['energy_max']):
		sys.stdout.write('----- {} -----\n'.format(sys._getframe().f_code.co_name))

		mask, message, suffix = self.get_energy_mask(
			energy_min=energy_min,energy_max=energy_max)

		lc = LightCurve(
			np.array(self.df[mask]['unixtime']),
			float(self.unixtime_offset),
			tbin=tbin,tstart=tstart,tstop=tstop)

		self.par_curve_stat = lc.set_curve_statistics(
			threshold_sigma=threshold_sigma,
			outpdf='%s/%s_cnthist_gaussfit.pdf' % (self.outdir,self.basename),
			title='%s Burst search (%s) %d-sigma' % (self.basename, message, threshold_sigma)
			)

		threshold_count = self.par_curve_stat['peak'] + threshold_sigma * self.par_curve_stat['sigma']
		self.par_curve_stat['threshold_sigma'] = threshold_sigma	
		self.par_curve_stat['threshold_count'] = threshold_count

		burst_mask = (lc.hist.y >= threshold_count)
		self.par_curve_stat['burst_time'] = lc.hist.x[burst_mask]

		print(self.par_curve_stat)

		lc.plot_burst_mask(
			burst_mask=burst_mask,
			outpdf='%s/%s_lc_%s_bst.pdf' % (self.outdir,self.basename,suffix),
			title='%s burst search (%s) %d-sigma' % (self.basename, message, threshold_sigma))	

		gti_start_index = np.argwhere(burst_mask[:-1] < burst_mask[1:]).squeeze() # [False,True] transition
		gti_stop_index = np.argwhere(burst_mask[:-1] > burst_mask[1:]).squeeze() # [False,True] transition
		print("gti_start_index",gti_start_index)
		print("gti_stop_index",gti_stop_index)
		print(lc.hist.x[gti_start_index])
		print(lc.hist.x[gti_start_index]-tbin*0.5)

	"""
	def fit_burst_curve(self,
			tbin=DICT_FIT_BURST_CURVE['tbin'],
			tstart=DICT_FIT_BURST_CURVE['tstart'],
			tstop=DICT_FIT_BURST_CURVE['tstop'],
			energy_min=DICT_FIT_BURST_CURVE['energy_min'],
			energy_max=DICT_FIT_BURST_CURVE['energy_max'],
			fit_nsigma=DICT_FIT_BURST_CURVE['fit_nsigma']):

		sys.stdout.write('----- {} -----\n'.format(sys._getframe().f_code.co_name))

		mask, message, suffix = self.get_energy_mask(energy_min=energy_min,energy_max=energy_max)

		lc = LightCurve(np.array(self.df[mask]['unixtime']),float(self.unixtime_offset),
			tbin=tbin,tstart=tstart,tstop=tstop)
		title = '%s (%s)' % (self.basename, message)

#		dict_initpar = {'peak':300,'sigma':10,'area':100,'c0':10,'c1':0}
		dict_initpar = {'peak':300,'sigma':10,'area':600,'c0':14,'c1':0,'fit_xmin':200,'fit_xmax':500}
		par = lc.fit_gauss_linear(dict_initpar,flag_error=True,fit_nsigma=fit_nsigma)
		print("---------")
		print(par)
		print("---------")

		model_x = lc.hist.x
		model_y = np.array([model_gauss_linear(x,peak=par['peak'],sigma=par['sigma'],area=par['area'],
			c0=par['c0'],c1=par['c1']) for x in model_x])	

		legend_text  = 'peak=%.1f+/-%.1f ch\n' % (par['peak'],par['peak_err'])
		legend_text += 'sigma=%.1f+/-%.1f ch\n' % (par['sigma'],par['sigma_err'])		
		legend_text += 'area=%.1f+/-%.1f counts\n' % (par['area'],par['area_err'])
		plot_fit_residual(
			lc.hist.x, lc.hist.y,
			model_x, model_y,			
			outpdf='%s/%s_bst_gaussfit.pdf' % (self.outdir,self.basename),
			hist_xerr=lc.hist.xerr,
			hist_yerr=lc.hist.yerr,
			xlabel='Time (sec) since %s' % datetime.fromtimestamp(self.unixtime_offset),
			ylabel='Counts / (%d sec)' % tbin,
			title='%s (%s)' % (self.basename, message),
			xlim=[par['fit_xmin'],par['fit_xmax']],
			axvline_values=[par['peak']-3*par['sigma'],par['peak']+3*par['sigma']],
			axvline_legends=[r"$-3\sigma$",r"$+3\sigma$"],			
			legend_text=legend_text)
		return par
	"""

	def get_burst_duration(self,par,tbin=1.0,tstart=0.0,tstop=500.0,
		energy_min=3.0,energy_max=8.0,
		linear_tbin_normalization=1.0):

		mask, message, suffix = self.get_energy_mask(energy_min=energy_min,energy_max=energy_max)

		unixtime_series = np.array(self.df[mask]['unixtime']) - float(self.unixtime_offset)
		lc = LightCurve(np.array(self.df[mask]['unixtime']),float(self.unixtime_offset),
			tbin=tbin,tstart=tstart,tstop=tstop)		
		print(unixtime_series)

		outpdf = '%s/%s_accumy.pdf' % (self.outdir,self.basename)
		lc.plot_accumulation(par,linear_tbin_normalization=linear_tbin_normalization,outpdf=outpdf)
		#title = '%s (%s)' % (self.basename, message)
		#lc.plot(outpdf='lc1.pdf',title=title,xlim=[0.0,500.0])		

	def write_to_fitsfile(self,output_fitsfile=None,config_file=None,overwrite=True):
		"""
		https://docs.astropy.org/en/stable/io/fits/usage/table.html
		"""
		sys.stdout.write('----- {} -----\n'.format(sys._getframe().f_code.co_name))

		if output_fitsfile == None:
			output_fitsfile = "%s/%s.evt" % (self.outdir,self.basename)
		elif os.path.exists(output_fitsfile):
			if overwrite:
				cmd = 'rm -f %s' % output_fitsfile
				print(cmd);os.system(cmd)
			else:
				raise FileExistsError("{} has alaredy existed.".format(output_fitsfile))

#		self.set_time_series()
#		print(type(self.time_series_utc.unix))
#		print(type(self.df['decisec']))

		column_minute = fits.Column(name='minute',format='B', unit='minute', array=self.df['minute'])
		column_sec = fits.Column(name='sec',format='B', unit='sec', array=self.df['sec'])
		column_decisec = fits.Column(name='decisec',format='I', unit='100 microsec', array=self.df['decisec'])		
		column_pha = fits.Column(name='pha',format='I', unit='channel', array=self.df['pha'])
		column_unixtime = fits.Column(name='unixtime',format='D', unit='sec', array=self.df['unixtime'])
		column_energy = fits.Column(name='energy',format='D', unit='MeV', array=self.df['energy_mev'])		

		column_defs = fits.ColDefs([column_minute,column_sec,column_decisec,column_pha,column_unixtime,column_energy])
		hdu = fits.BinTableHDU.from_columns(column_defs,name='EVENTS')

		dict_keywords = {
			'DET_ID':[self.detid_str,'Detector_ID'],
			'YYYYMMDD':[self.yyyymmdd_jst,'Year, month, and day in JST of the file'],			
			'Hour':[self.hour_jst,'Hour in JST of the file']
			}
		for keyword in dict_keywords.keys():
			hdu.header[keyword] = dict_keywords[keyword][0]
			hdu.header.comments[keyword] = dict_keywords[keyword][1]

		if config_file != None:
			self.set_config_file(config_file)
			for keyword in self.config.dict_keywords.keys():
				hdu.header[keyword] = self.config.dict_keywords[keyword]

		hdu.header['comment'] = 'unixtime is UTC, while minute, sec, decisec columns and the file name are JST.'
		hdu.header['history'] = 'created at {} JST'.format(Time.now().to_datetime(tz_tokyo))
		hdu.writeto(output_fitsfile)

	def write_to_yamlfile(self):
		yamlfile = '%s/%s_process.yaml' % (self.outdir,self.basename)
		with open(yamlfile, "w") as wf:
		    yaml.dump(self.param, wf,default_flow_style=False)		

	def __str__(self):
		dump  = 'filepath: {}\n'.format(self.filepath)
		dump += 'filetype: {}\n'.format(self.filetype)		
		dump += '\n'
		return dump

	def show_summary(self,show_head_nrows=8):
		sys.stdout.write('----- {} -----\n'.format(sys._getframe().f_code.co_name))		
		sys.stdout.write(str(self))
		pd.options.display.float_format = '{:.5f}'.format
		pd.options.display.precision = 5
		print(self.df.head(show_head_nrows))

class PhaSpectrum():
	"""
	Spectrum
	"""
	def __init__(self, pha_series, nbins, pha_min, pha_max):
		self.nbins = nbins
		self.pha_min = pha_min
		self.pha_max = pha_max 

		self.hist = Hist1D(nbins=self.nbins,xlow=self.pha_min-0.5,xhigh=self.pha_max-0.5)
		self.hist.fill(pha_series)		

	def plot(self,outpdf,title='',xlim=[1,PHA_SPECTRUM_NBINS]):
		plot_histogram(
			self.hist.x,self.hist.y,
			outpdf=outpdf,hist_yerr=None,
			xlabel='ADC channel (pha)',ylabel='Counts / bin',title=title,
			flag_xlog=False,flag_ylog=True,xlim=xlim)		

	def fit_gauss_linear(self,par_init,flag_error=True,fit_nsigma=3):
		par = fit_gauss_linear(self.hist.x,self.hist.y,
			par_init,error=self.hist.yerr,fit_nsigma=fit_nsigma)
		return par 

class EnergySpectrum():
	"""
	Energy Spectrum
	"""
	def __init__(self, energy_series, nbins, energy_min, energy_max):
		self.nbins = nbins
		self.energy_min = energy_min
		self.energy_max = energy_max 

		self.hist = Hist1D(nbins=self.nbins,xlow=self.energy_min,xhigh=self.energy_max)
		self.hist.fill(energy_series)		

	def plot(self,outpdf,title='',
		xlim=[ENERGY_SPECTRUM_MIN,ENERGY_SPECTRUM_MAX]):
		plot_histogram(
			self.hist.x,self.hist.y,
			outpdf=outpdf,hist_yerr=None,
			xlabel='Energy (MeV)',ylabel='Counts / bin',title=title,
			flag_xlog=False,flag_ylog=True,xlim=xlim)		

class LightCurve():
	"""
	Light curve
	"""
	def __init__(self,unixtime_series,unixtime_offset,tbin=8.0,tstart=0.0,tstop=3600.0):
		self.unixtime_series = unixtime_series
		self.unixtime_offset = unixtime_offset 
		self.tbin = tbin
		self.tstart = tstart
		self.tstop = tstop 

		self.nbins = round((self.tstop - self.tstart)/self.tbin)
		self.hist = Hist1D(nbins=self.nbins,xlow=self.tstart,xhigh=self.tstop)
		self.hist.fill(self.unixtime_series-self.unixtime_offset)

		self.time_offset_str = datetime.fromtimestamp(self.unixtime_offset)

	def plot(self,outpdf,title='',xlim=[0.0,3600.0]):

		plot_histogram(
			self.hist.x,self.hist.y,
			outpdf=outpdf,hist_yerr=self.hist.yerr,
			xlabel='Time (sec) since %s JST' % self.time_offset_str,
			ylabel='Counts / (%d sec)' % self.tbin,title=title,
			flag_xlog=False,flag_ylog=False,xlim=xlim)	

	def plot_burst_mask(self,burst_mask,outpdf,title='',xlim=[0.0,3600.0]):

		plot_histogram(
			self.hist.x,self.hist.y,
			outpdf=outpdf,hist_yerr=self.hist.yerr,
			xlabel='Time (sec) since %s JST' % self.time_offset_str,
			ylabel='Counts / (%d sec)' % self.tbin,title=title,
			flag_xlog=False,flag_ylog=False,xlim=xlim,
			burst_mask=burst_mask)	

	def set_curve_statistics(self,outpdf,title='',threshold_sigma=5.0):
		sys.stdout.write('----- {} -----\n'.format(sys._getframe().f_code.co_name))

		cnt_xlim_max = round(1.2*max(self.hist.y))
		self.cnt_hist = Hist1D(nbins=int(cnt_xlim_max),
			xlow=-0.5,xhigh=cnt_xlim_max-0.5)
		self.cnt_hist.fill(self.hist.y)

		peak = np.mean(self.hist.y)
		sigma = np.std(self.hist.y)
		area = np.sum(self.hist.y)
		dict_initpar = {'peak':peak,'sigma':sigma,'area':area,
			'fit_xmin':-0.5,'fit_xmax':cnt_xlim_max}
		par = fit_gauss(self.cnt_hist.x,self.cnt_hist.y,
			par_init=dict_initpar,error=None,fit_nsigma=5)
		par['tbin'] = self.tbin
		par['tstart'] = self.tstart
		par['tstop'] = self.tstop

		model_x = self.cnt_hist.x
		model_y = np.array([model_gauss(x,peak=par['peak'],sigma=par['sigma'],area=par['area']) for x in model_x])			

		legend_text  = 'peak=%.2f+/-%.2f\n' % (par['peak'],par['peak_err'])
		legend_text += 'sigma=%.2f+/-%.2f\n' % (par['sigma'],par['sigma_err'])	
		legend_text += 'area=%.1f+/-%.1f\n' % (par['area'],par['area_err'])
		plot_fit_residual(
			self.cnt_hist.x, self.cnt_hist.y,
			model_x, model_y,			
			hist_xerr=self.cnt_hist.xerr,
			hist_yerr=self.cnt_hist.yerr,
			outpdf=outpdf, 
			xlabel='Counts / (%.1f sec bin)' % self.tbin,
			ylabel='Number of bins',
			title=title, 
			xlim=[-0.5,cnt_xlim_max-0.5],
			axvline_values=[par['peak']-threshold_sigma*par['sigma'],par['peak']+threshold_sigma*par['sigma']],
			axvline_legends=[r"$-%d\sigma$" % threshold_sigma,r"$+%d\sigma$ (prob.=%.2e)" % (threshold_sigma,get_norm_probability(threshold_sigma))],
			legend_text=legend_text,legend_loc='upper right')
		return par

	def fit_gauss_linear(self,par_init,flag_error=True,fit_nsigma=3):
		par = fit_gauss_linear(self.hist.x,self.hist.y,par_init,error=self.hist.yerr,fit_nsigma=fit_nsigma)
		return par 

	def plot_accumulation(self,par,linear_tbin_normalization=1.0,outpdf='accum.pdf'):

		model_x = self.hist.x
		model_y = np.array([model_linear(x,c0=par['c0'],c1=par['c1'])/linear_tbin_normalization for x in model_x])
		#print(model_x)
		#print(model_y)

		self.hist.accum_y = (self.hist.y - model_y).cumsum()
		#print(self.hist.accum_y)

		ratio = self.hist.accum_y/self.hist.accum_y[-1]
		at10percent = self.hist.x[ratio >= 0.1][0]
		at90percent = self.hist.x[ratio >= 0.9][0]
		print(at10percent,at90percent)
		print(at90percent-at10percent)

		title = 'accumulation'
		plot_histogram(
			self.hist.x,self.hist.accum_y,
			outpdf=outpdf,hist_yerr=None,
			xlabel='Time (sec) since %s JST' % self.time_offset_str,
			ylabel='Number of events',title=title,
			flag_xlog=False,flag_ylog=False,xlim=None)	



