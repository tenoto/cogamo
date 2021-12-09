# -*- coding: utf-8 -*-

import os 
import re
import sys
import yaml
import numpy as np 
import pandas as pd 

from astropy.time import Time
from astropy.io import fits

from datetime import datetime, timedelta, timezone
tz_tokyo = timezone(timedelta(hours=+9), 'Asia/Tokyo')
tz_utc = timezone(timedelta(hours=0), 'UTC')

import matplotlib.pylab as plt 
import matplotlib.gridspec as gridspec

from iminuit import Minuit
from probfit import Chi2Regression

##################################################
# Default value 
##################################################

DEFAULT_PHA_SPECTRUM_NBINS = 2**10

DEFAULT_ENERGY_SPECTRUM_NBINS = 2**10
DEFAULT_ENERGY_SPECTRUM_MIN = 0.0
DEFAULT_ENERGY_SPECTRUM_MAX = 12.0

K40_ENERGY_MEV = 1.46083
Tl208_ENERGY_MEV = 2.61453

dict_par_init_Tl208 = {'peak':236,'sigma':7,'area':2651,'c0':798.0,'c1':-3,'MeV':Tl208_ENERGY_MEV,'pha_min':200,'pha_max':264,'nbins':64,'xlim':[200,264],'name':'Tl208'}
dict_par_init_K40 = {'peak':132,'sigma':5,'area':18025,'c0':3731,'c1':-21.0,'MeV':K40_ENERGY_MEV,'pha_min':100,'pha_max':164,'nbins':64,'xlim':[100,164],'name':'K40'}

##################################################
# General function and classes 
##################################################

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
		flag_xlog=False,flag_ylog=False,xlim=None):
	fig, ax = plt.subplots(1,1, figsize=(11.69,8.27)) # A4 size, inich unit 
	fontsize = 18 	

	if hist_yerr is not None:
		plt.errorbar(hist_x,hist_y,marker='',drawstyle='steps-mid')		
		plt.errorbar(hist_x,hist_y,yerr=hist_yerr,marker='',drawstyle='steps-mid')
	else:
		plt.errorbar(hist_x,hist_y,marker='',drawstyle='steps-mid')

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

def plot_fit_residual(hist_x,hist_y,model_x,model_y,outpdf,
		hist_xerr=None,hist_yerr=None,
		xlabel='X title',ylabel='Y title',title='Title',
		flag_xlog=False,flag_ylog=False,xlim=None,
		fit_xmin=None,fit_xmax=None,legend_text=''):

	fig, axs = plt.subplots(2,1,figsize=(11.69,8.27), # A4 size, inich unit 
		sharex=True,gridspec_kw={'hspace':0},tight_layout=True)
	gs = gridspec.GridSpec(2,1,height_ratios=[3,1])
	gs.update(hspace=0.0)

	fontsize = 18
	plt.tight_layout()
	plt.tick_params(labelsize=fontsize)
	plt.rcParams["font.family"] = "serif"
	plt.rcParams["mathtext.fontset"] = "dejavuserif"	

	axs[0] = plt.subplot(gs[0,0])
	axs[0].errorbar(hist_x,hist_y,xerr=hist_xerr,yerr=hist_yerr,
		marker='o',linestyle='None',color='k')
	axs[0].plot(model_x,model_y,c='red',drawstyle='steps-mid')
	axs[0].set_ylabel(ylabel, fontsize=fontsize)
	axs[0].set_xlim(xlim)
	axs[0].set_title(title, fontsize=fontsize)	
	axs[0].get_xaxis().set_visible(False)
	axs[0].legend(title=legend_text,loc='lower left',fontsize=fontsize)

	axs[1] = plt.subplot(gs[1])
	axs[1].errorbar(hist_x,(hist_y-model_y)/model_y,yerr=hist_yerr/model_y,
		marker='o',linestyle='None',color='k')
	axs[1].axhline(y=0,ls='--',color='r')
	axs[1].set_ylabel('(data-model)/model', fontsize=fontsize)
	axs[1].set_xlabel(xlabel, fontsize=fontsize)
	axs[1].set_xlim(xlim)

	for ax in axs:
		ax.minorticks_on()
		#ax.grid(True)
		#ax.grid(axis='both',which='major', linestyle='--', color='#000000')
		#ax.grid(axis='both',which='minor', linestyle='--')	
		ax.tick_params(axis="both", which='major', direction='in', length=5,labelsize=fontsize)
		ax.tick_params(axis="both", which='minor', direction='in', length=3,labelsize=fontsize)
		ax.axvline(fit_xmin,ls='--')
		ax.axvline(fit_xmax,ls='--')

	plt.savefig(outpdf)	

def model_gauss_linear(x, peak, sigma, area, c0=0.0, c1=0.0):
    return area * np.exp(-0.5*(x-peak)**2/sigma**2)/(np.sqrt(2*np.pi)*sigma) + c0 + c1 * x

def model_linear(x, c0=0.0, c1=0.0):
    return c0 + c1 * x

def fit_gauss_linear(x,y,par_init,error=None,fit_nsigma=3):

	if error is not None:
		chi2reg = Chi2Regression(model_gauss_linear,x,y,error=error)
	else:
		chi2reg = Chi2Regression(model_gauss_linear,x,y)	

	fit = Minuit(chi2reg, 
		peak=par_init['peak'],sigma=par_init['sigma'],area=par_init['area'], 
		c0=par_init['c0'], c1=par_init['c1'])
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

class Pipeline():
	def __init__(self):
		sys.stdout.write("cogamo pipeline\n")

	def process_eventdata(self,filepath):

		# initial setup 
		param = {}

		# read raw csv file 
		evt = EventData(filepath)
		if evt.filetype != 'rawcsv':
			sys.stdout.write('Input file is not the rawcsv file format.')
			return 0		

		# prepare outputfile 
		outdir = 'out/%s' % evt.basename
		cmd = 'mkdir -p %s' % outdir 
		os.system(cmd)

		# get energy calibration informatio: two line peak channel 
		evt.extract_pha_spectrum(outdir=outdir)
		par_Tl208 = evt.fit_pha_spectrum_lines(dict_par_init_Tl208,outdir=outdir)
		par_K40 = evt.fit_pha_spectrum_lines(dict_par_init_K40,outdir=outdir)

		for par in [par_K40,par_Tl208]:
			for strtype in ['','_err']:
				param['%s_peak%s' % (par['name'],strtype)] = par['peak%s' % strtype]			
				param['%s_sigma%s' % (par['name'],strtype)] = par['sigma%s' % strtype]
				param['%s_area%s' % (par['name'],strtype)] = par['area%s' % strtype]
				param['%s_c0%s' % (par['name'],strtype)] = par['c0%s' % strtype]			
				param['%s_c1%s' % (par['name'],strtype)] = par['c1%s' % strtype]						
			param['%s_name' % par['name']] = par['name']						
			param['%s_MeV' % par['name']] = par['MeV']									

		outpdf = '%s/%s_energycal.pdf' % (outdir,evt.basename)
		pha2mev_c0, pha2mev_c1 = self.get_pha2mev_param(
			mev_array=np.array([K40_ENERGY_MEV,Tl208_ENERGY_MEV]),
			pha_array=np.array([param['K40_peak'],param['Tl208_peak']]),
			pha_err_array=np.array([param['K40_peak_err'],param['Tl208_peak_err']]),
			title=evt.basename,outpdf=outpdf)
		param['pha2mev_c0'] = pha2mev_c0
		param['pha2mev_c1'] = pha2mev_c1

		# energy calibration 
		evt.set_energy_series(pha2mev_c0=pha2mev_c0,pha2mev_c1=pha2mev_c1)

		# time assignment 
		evt.set_time_series()

		# energy spectrum 
		evt.extract_energy_spectrum(outdir=outdir)

		evt.extract_curve(tbin=5.0,energy_min=3.0,energy_max=None,outdir=outdir)

		# if burst was detected
		#par = evt.fit_burst_curve(tbin=5.0,tstart=0.0,tstop=500.0,energy_min=3.0,energy_max=None)
		par = evt.fit_burst_curve(tbin=5.0,tstart=0.0,tstop=500.0,energy_min=3.0,energy_max=None,outdir=outdir)		

		evt.get_burst_duration(par,tbin=1.0,tstart=par['fit_xmin'],tstop=par['fit_xmax'],
			energy_min=3.0,energy_max=None,linear_tbin_normalization=5.0,outdir=outdir)

		# output 
		output_fitsfile = '%s/%s_proc.evt' % (outdir,evt.basename)
		evt.write_to_fitsfile(output_fitsfile=output_fitsfile,config_file=None)

		yamlfile = '%s/%s_process.yaml' % (outdir,evt.basename)
		with open(yamlfile, "w") as wf:
		    yaml.dump(param, wf,default_flow_style=False)	

	def get_pha2mev_param(self,mev_array,pha_array,pha_err_array,title,outpdf):

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

		self.set_filetype()
		self.open_file()

	def set_filetype(self):
		if re.fullmatch(r'\d{3}_\d{8}_\d{2}.csv', self.filename):
			self.filetype = 'rawcsv'
			self.detid_str, self.yyyymmdd_jst, self.hour_jst = self.basename.split("_")		
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

	def set_time_series(self):		
		"""
		the standard unix time does not have enough accuracy below 1 second.
		Sub-second time stamp is handled by another column
		"""
		sys.stdout.write('----- {} -----\n'.format(sys._getframe().f_code.co_name))

		year = self.yyyymmdd_jst[0:4]
		month = self.yyyymmdd_jst[4:6]
		day = self.yyyymmdd_jst[6:8]		
		str_time = '%04d-%02d-%02dT%02d:' % (int(year),int(month),int(day),int(self.hour_jst))
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
			nbins=DEFAULT_PHA_SPECTRUM_NBINS,pha_min=0.0,pha_max=DEFAULT_PHA_SPECTRUM_NBINS,
			xlim=[1,DEFAULT_PHA_SPECTRUM_NBINS],outdir='./'):
		"""Extract spectrum of pha
		nbins: pha bin size
		"""
		sys.stdout.write('----- {} -----\n'.format(sys._getframe().f_code.co_name))

		phaspec = PhaSpectrum(self.df['pha'],nbins=nbins,pha_min=pha_min,pha_max=pha_max)
		outpdf = '%s/%s_phaspec.pdf' % (outdir,self.basename)
		phaspec.plot(outpdf=outpdf,title=self.basename,xlim=xlim)
		return phaspec

	def extract_energy_spectrum(self,
			nbins=DEFAULT_ENERGY_SPECTRUM_NBINS,
			energy_min=DEFAULT_ENERGY_SPECTRUM_MIN,
			energy_max=DEFAULT_ENERGY_SPECTRUM_MAX,
			xlim=[DEFAULT_ENERGY_SPECTRUM_MIN,DEFAULT_ENERGY_SPECTRUM_MAX],
			outdir='./'):
		"""Extract spectrum of pha
		nbins: pha bin size
		"""
		sys.stdout.write('----- {} -----\n'.format(sys._getframe().f_code.co_name))

		enespec = EnergySpectrum(self.df['energy_mev'],
			nbins=nbins,energy_min=energy_min,energy_max=energy_max)
		outpdf = '%s/%s_energyspec.pdf' % (outdir,self.basename)
		enespec.plot(outpdf=outpdf,title=self.basename,xlim=xlim)
		return enespec

	def fit_pha_spectrum_lines(self,dict_par_init,outdir='./'):
		"""Extract spectrum of pha
		nbins: pha bin size
		"""
		sys.stdout.write('----- {} -----\n'.format(sys._getframe().f_code.co_name))

		phaspec = PhaSpectrum(self.df['pha'],nbins=dict_par_init['nbins'],
			pha_min=dict_par_init['pha_min'],pha_max=dict_par_init['pha_max'])
		outpdf = '%s/%s_%s.pdf' % (outdir,self.basename,dict_par_init['name'])

		par = phaspec.fit_gauss_linear(dict_par_init,flag_error=True)
		par['name'] = dict_par_init['name']
		par['MeV'] = dict_par_init['MeV']		

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
			xlim=[dict_par_init['pha_min'],dict_par_init['pha_max']],
			fit_xmin=par['fit_xmin'],fit_xmax=par['fit_xmax'],
			legend_text=legend_text)
		return par


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

	def extract_curve(self,tbin=1.0,tstart=0.0,tstop=3600.0,energy_min=3.0,energy_max=8.0,outdir='./'):
		""" Energy in MeV unit
		"""
		mask, message, suffix = self.get_energy_mask(energy_min=energy_min,energy_max=energy_max)

		lc = LightCurve(np.array(self.df[mask]['unixtime']),float(self.unixtime_offset),tbin=tbin,tstart=tstart,tstop=tstop)
		title = '%s (%s)' % (self.basename, message)
		outpdf = '%s/%s_lc_%s.pdf' % (outdir,self.basename,suffix)
		lc.plot(outpdf=outpdf,title=title,xlim=[0.0,500.0])		

	def fit_burst_curve(self,tbin=8.0,tstart=0.0,tstop=3600.0,energy_min=3.0,energy_max=8.0,outdir='./'):
		""" Energy in MeV unit
		"""
		mask, message, suffix = self.get_energy_mask(energy_min=energy_min,energy_max=energy_max)

		lc = LightCurve(np.array(self.df[mask]['unixtime']),float(self.unixtime_offset),
			tbin=tbin,tstart=tstart,tstop=tstop)
		title = '%s (%s)' % (self.basename, message)

#		dict_par_init = {'peak':300,'sigma':10,'area':100,'c0':10,'c1':0}
		dict_par_init = {'peak':300,'sigma':10,'area':600,'c0':14,'c1':0,'fit_xmin':200,'fit_xmax':500}
		#dict_par_init = {'peak':300,'sigma':10,'area':600,'c0':3,'c1':0,'fit_xmin':200,'fit_xmax':500}
		par = lc.fit_gauss_linear(dict_par_init,flag_error=True,fit_nsigma=6)
		print("---------")
		print(par)
		print("---------")

		model_x = lc.hist.x
		model_y = np.array([model_gauss_linear(x,peak=par['peak'],sigma=par['sigma'],area=par['area'],
			c0=par['c0'],c1=par['c1']) for x in model_x])	

		"""
		legend_text = '%s (fit range: %d sigma)\n' % (self.basename,par['fit_nsigma'])
		legend_text += 'peak=%.1f+/-%.1f ch\n' % (par['peak'],par['peak_err'])
		legend_text += 'sigma=%.1f+/-%.1f ch\n' % (par['sigma'],par['sigma_err'])		
		legend_text += 'area=%.1f+/-%.1f counts\n' % (par['area'],par['area_err'])
		legend_text += 'c0=%.1f+/-%.1f counts\n' % (par['c0'],par['c0_err'])		
		legend_text += 'c1=%.1f+/-%.1f counts\n' % (par['c1'],par['c1_err'])			
		legend_text += 'resolution=%.1f %%' % (100.0*get_energy_resolution(par['peak'],par['sigma']))
		"""
		legend_text = ""
		outpdf = '%s/%s_bst_gaussfit.pdf' % (outdir,self.basename)
		plot_fit_residual(
			lc.hist.x,lc.hist.y,
			model_x,model_y,			
			outpdf=outpdf,
			hist_xerr=lc.hist.xerr,hist_yerr=lc.hist.yerr,
			xlabel='Channel',ylabel='Counts / (%d sec)' % tbin,
			title=self.basename,
			xlim=[tstart,tstop],
			fit_xmin=par['fit_xmin'],fit_xmax=par['fit_xmax'],
			legend_text=legend_text)
		return par

	def get_burst_duration(self,par,tbin=1.0,tstart=0.0,tstop=500.0,energy_min=3.0,energy_max=8.0,
		linear_tbin_normalization=1.0,outdir='./'):

		mask, message, suffix = self.get_energy_mask(energy_min=energy_min,energy_max=energy_max)

		unixtime_series = np.array(self.df[mask]['unixtime']) - float(self.unixtime_offset)
		lc = LightCurve(np.array(self.df[mask]['unixtime']),float(self.unixtime_offset),
			tbin=tbin,tstart=tstart,tstop=tstop)		
		print(unixtime_series)

		outpdf = '%s/%s_accumy.pdf' % (outdir,self.basename)
		lc.plot_accumulation(par,linear_tbin_normalization=linear_tbin_normalization,outpdf=outpdf)
		#title = '%s (%s)' % (self.basename, message)
		#lc.plot(outpdf='lc1.pdf',title=title,xlim=[0.0,500.0])		

	def write_to_fitsfile(self,output_fitsfile=None,config_file=None,overwrite=True):
		"""
		https://docs.astropy.org/en/stable/io/fits/usage/table.html
		"""
		sys.stdout.write('----- {} -----\n'.format(sys._getframe().f_code.co_name))

		if output_fitsfile == None:
			output_fitsfile = "{}.evt".format(self.basename)
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

	def plot(self,outpdf,title='',xlim=[1,DEFAULT_PHA_SPECTRUM_NBINS]):
		plot_histogram(
			self.hist.x,self.hist.y,
			outpdf=outpdf,hist_yerr=None,
			xlabel='ADC channel (pha)',ylabel='Counts / bin',title=title,
			flag_xlog=False,flag_ylog=True,xlim=xlim)		

	def fit_gauss_linear(self,par_init,flag_error=True,fit_nsigma=3):
		par = fit_gauss_linear(self.hist.x,self.hist.y,par_init,error=self.hist.yerr,fit_nsigma=fit_nsigma)
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
		xlim=[DEFAULT_ENERGY_SPECTRUM_MIN,DEFAULT_ENERGY_SPECTRUM_MAX]):
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
			xlabel='Time since %s (sec)' % self.time_offset_str,
			ylabel='Counts / (%d sec)' % self.tbin,title=title,
			flag_xlog=False,flag_ylog=False,xlim=xlim)	

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
			xlabel='Time since %s (sec)' % self.time_offset_str,
			ylabel='Number of events',title=title,
			flag_xlog=False,flag_ylog=False,xlim=None)	



