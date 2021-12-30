# -*- coding: utf-8 -*-
# Last modified 2021-12-27

import os 
import re
import sys
import yaml
import glob 
import numpy as np 
import pandas as pd 
from scipy.stats import norm

from astropy.time import Time
from astropy.io import fits

from datetime import datetime, timedelta, timezone
tz_tokyo = timezone(timedelta(hours=+9), 'Asia/Tokyo')
tz_utc = timezone(timedelta(hours=0), 'UTC')

import matplotlib.pylab as plt 
import matplotlib.gridspec as gridspec

from iminuit import Minuit
from probfit import Chi2Regression

import PyPDF2

MAX_PHA_CHANNEL = 2**10 - 1 		

"""
Response

  No. Type     EXTNAME      BITPIX Dimensions(columns)      PCOUNT  GCOUNT
 
   0  PRIMARY                 16     0                           0    1
   1  BINTABLE MATRIX          8     8206(6) 592                 0    1
 
      Column Name                Format     Dims       Units     TLMIN  TLMAX
      1 ENERG_LO                   1E
      2 ENERG_HI                   1E
      3 N_GRP                      1I
      4 F_CHAN                     1I
      5 N_CHAN                     1I
      6 MATRIX                     2048E
 
   2  BINTABLE EBOUNDS         8     10(3) 2048                  0    1
 
      Column Name                Format     Dims       Units     TLMIN  TLMAX
      1 CHANNEL                    1I
      2 E_MIN                      1E
      3 E_MAX                      1E
 
[MATRIX]
Channel ENERG_LO ENERG_HI
(ch)    (keV)    (keV)
1       40       45 
2      	45		 50
3		50		 55
...
592		40800	 41000

[EBOUNDS]
Channel ENERG_LO ENERG_HI
(ch)    (keV)    (keV)
0       40       60
1      	60		 80
2		80		 100
...
2046	40960	 40980
2047	40980	 41000

ENERG_LO = 40 keV + 20 keV * i 
ENERG_HI = 60 keV + 20 keV * i 
"""

#DICT_INITPAR_K40 = {'name':'K40','MeV':1.46083,'peak':132,'sigma':5,'area':18025,'c0':3731,'c1':-21.0,'pha_min':100,'pha_max':164,'binning':2,'xlim':[100,164]}
DICT_INITPAR_K40 = {'name':'K40','MeV':1.46083,'peak':132,'sigma':5,'area':18025,'c0':3731,'c1':-21.0,'pha_min':100,'pha_max':184,'binning':2,'xlim':[100,184]}
#DICT_INITPAR_TL208 = {'name':'Tl208','MeV':2.61453,'peak':236,'sigma':7,'area':2651,'c0':798.0,'c1':-3,'pha_min':190,'pha_max':284,'binning':2,'xlim':[190,284],}
#DICT_INITPAR_TL208 = {'name':'Tl208','MeV':2.61453,'peak':250,'sigma':7,'area':2651,'c0':798.0,'c1':-3,'pha_min':210,'pha_max':300,'binning':2,'xlim':[210,300],}
DICT_INITPAR_TL208 = {'name':'Tl208','MeV':2.61453,'peak':220,'sigma':7,'area':2651,'c0':798.0,'c1':-3,'pha_min':190,'pha_max':284,'binning':2,'xlim':[190,284],}
GAMMA_LINES = [DICT_INITPAR_K40,DICT_INITPAR_TL208]

def model_gauss(x, peak, sigma, area):
    return area * np.exp(-0.5*(x-peak)**2/sigma**2)/(np.sqrt(2*np.pi)*sigma) 

def model_gauss_linear(x, peak, sigma, area, c0=0.0, c1=0.0):
    return area * np.exp(-0.5*(x-peak)**2/sigma**2)/(np.sqrt(2*np.pi)*sigma) + c0 + c1 * x

def model_linear(x, c0=0.0, c1=0.0):
    return c0 + c1 * x

def get_energy_resolution(peak,sigma):
	# assuming the offset is zero
	return 2.35*sigma/peak 

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

def extract_xspec_pha(energy_keV_array,outpha,exposure):
	number_of_channel = 2048
	energy_min = 40 # keV 
	energy_step = 20 # keV 

	energy_max = energy_min + energy_step * number_of_channel
	hist = Hist1D(nbins=number_of_channel,xlow=energy_min,xhigh=energy_max)
	hist.fill(energy_keV_array)

	col_channel = fits.Column(name='CHANNEL',format="J", array=np.array([i for i in range(0,number_of_channel)]))
	col_counts = fits.Column(name='COUNTS',format="J",unit="count",	array=hist.y)
	hdu1 = fits.BinTableHDU.from_columns([col_channel,col_counts])
	hdu1.name = 'SPECTRUM'

	prhdu = fits.PrimaryHDU()
	hdu = fits.HDUList([prhdu, hdu1])

	hdu1.header['HDUCLASS'] = 'OGIP'
	hdu1.header['HDUCLAS1'] = 'SPECTRUM'
	hdu1.header['HDUVERS1'] = '1.2.0   '   	
	hdu1.header['HDUVERS'] = '1.2.0   '
	hdu1.header['TLMIN1']  = 0   
	hdu1.header['TLMAX1']  = 2047
	hdu1.header['TELESCOP']= 'thdr'
	hdu1.header['INSTRUME']= 'CogamoCsI'
	hdu1.header['FILTER']  = 'UNKNOWN '
	hdu1.header['EXPOSURE'] = exposure 
	hdu1.header['AREASCAL'] = 1.0 
	hdu1.header['BACKFILE'] = 'NONE'
	hdu1.header['BACKSCAL'] = 1.0
	hdu1.header['ANCRFILE'] = 'NONE'
	hdu1.header['CORRFILE'] = 'NONE'
	hdu1.header['CORRSCAL'] = 1.0
	hdu1.header['PHAVERSN'] = '1992a'
	hdu1.header['CHANTYPE'] = 'PI'
	hdu1.header['POISSERR'] = True 
	hdu1.header['STAT_ERR'] = 0
	hdu1.header['SYS_ERR'] = 0	

	"""
	for i in range(0,2):
	  	for keyword,usrsetup in setup['keywords_common'].items():
			hdu[i].header[keyword] = usrsetup[0]
			hdu[i].header.comments[keyword] = usrsetup[1]
		hdu[i].header["FILENAME"] = output_rmffile
		hdu[i].header.comments[keyword] = "Filename"				
	"""		
	i = 1
	hdu[i].header["DETCHANS"] = number_of_channel
	#hdu[i].header.comments[keyword] = "total number of detector channels"

	hdu.writeto(outpha)  

class Hist1D(object):
	"""Represents 1D histogram. 
	"""
	def __init__(self, nbins, xlow, xhigh):
		self.nbins = nbins
		self.xlow  = xlow
		self.xhigh = xhigh
		self.y, edges = np.histogram([], 
			bins=nbins, range=(xlow, xhigh))
		self.x = (edges[:-1] + edges[1:]) / 2.
		self.xerr = (-edges[:-1] + edges[1:]) / 2.
		self.yerr = None

	def fill(self, arr):
		y, edges = np.histogram(arr, bins=self.nbins, range=(self.xlow, self.xhigh))
		self.y += y
		self.yerr = np.sqrt(self.y)

	@property
	def data(self):
		return self.x, self.y

	def plot(self,outpdf,
		xlabel='X title',ylabel='Y title',title='Title',		
		flax_xerr=False,flag_yerr=False,
		flag_xlog=False,flag_ylog=False,
		xlim=None,
		mask=None,
		axvline_values=None,axvline_legends=None,
		axhline_values=None,axhline_legends=None,
		legend_title="",legend_loc="upper right",legend_fontsize=16):

		fig, ax = plt.subplots(1,1,figsize=(11.69,8.27)) # A4 size, inich unit 
		fontsize = 18 	

		if flag_yerr:
			plt.errorbar(self.x,self.y,marker='',drawstyle='steps-mid',color='k')		
			plt.errorbar(self.x,self.y,yerr=self.yerr,marker='',drawstyle='steps-mid',color='k')
		else:
			plt.errorbar(self.x,self.y,marker='',drawstyle='steps-mid',color='k')

		if mask is not None:
			if flag_yerr is not None:
				plt.errorbar(self.x[mask],self.y[mask],marker='',drawstyle='steps-mid',color='r')		
				plt.errorbar(self.x[mask],self.y[mask],yerr=self.yerr[mask],marker='',drawstyle='steps-mid',color='r')
			else:
				plt.errorbar(self.x[mask],self.y[mask],marker='',drawstyle='steps-mid',color='r')		

		plt.xlabel(xlabel, fontsize=fontsize)
		plt.ylabel(ylabel, fontsize=fontsize)
		plt.title(title,fontsize=fontsize)
		if flag_xlog: plt.xscale('log')				
		if flag_ylog: plt.yscale('log')
		if xlim is not None: plt.xlim(xlim)
		
		ax.minorticks_on()
		ax.grid(True)
		ax.grid(axis='both',which='major', linestyle='--', color='#000000')
		ax.grid(axis='both',which='minor', linestyle='--')	
		ax.tick_params(axis="both", which='major', direction='in', length=5)
		ax.tick_params(axis="both", which='minor', direction='in', length=3)

		if axvline_values is not None:
			if axvline_legends is None:
				for value in axvline_values:
					ax.axvline(value,ls='--')				
			else:
				for i in range(len(axvline_values)):
					ax.axvline(axvline_values[i],ls='--',label=axvline_legends[i])

		if axhline_values is not None:
			if axhline_legends is None:
				for value in axhline_values:
					ax.axhline(value,ls='--')				
			else:
				for i in range(len(axhline_values)):
					ax.axhline(axhline_values[i],ls='--',label=axhline_legends[i])

		legend = ax.legend(title=legend_title,loc=legend_loc,fontsize=legend_fontsize)

		plt.tick_params(labelsize=fontsize)
		plt.rcParams["font.family"] = "serif"
		plt.rcParams["mathtext.fontset"] = "dejavuserif"	
		plt.tight_layout()

		plt.savefig(outpdf)			

	def fit_gauss(self,par_init,flag_error=False,fit_nsigma=3):

		mask = np.logical_and(self.x >= par_init["fit_xmin"], self.x <= par_init["fit_xmax"])

		if flag_error:
			chi2reg = Chi2Regression(model_gauss,
				self.x[mask],self.y[mask],error=self.yerr[mask])
		else:
			chi2reg = Chi2Regression(model_gauss,
				self.x[mask],self.y[mask])	

		fit = Minuit(chi2reg, 
			peak=par_init['peak'],sigma=par_init['sigma'],area=par_init['area'],
			limit_peak=(0,None),limit_sigma=(0,None),limit_area=(0,None))
		fit.migrad()
		fit.minos()
		fit.print_param()

		peak = fit.values[0]
		sigma = fit.values[1]
		fit_xmin = peak - fit_nsigma * sigma 
		fit_xmax = peak + fit_nsigma * sigma 

		mask = np.logical_and(self.x >= fit_xmin, self.x <= fit_xmax)
		if flag_error:
			chi2reg = Chi2Regression(model_gauss,
				self.x[mask],self.y[mask],error=self.yerr[mask])
		else:
			chi2reg = Chi2Regression(model_gauss,
				self.x[mask],self.y[mask])

		fit = Minuit(chi2reg, 
			peak=par_init['peak'],sigma=par_init['sigma'],area=par_init['area'],
			limit_peak=(0,None),limit_sigma=(0,None),limit_area=(0,None))
		fit.migrad()
		fit.minos() 
		fit.print_param()

		par = {
			'peak':fit.values[0],
			'sigma':fit.values[1],
			'area':fit.values[2],
			'peak_err':fit.errors[0],
			'sigma_err':fit.errors[1],
			'area_err':fit.errors[2],
			'fval':fit.fval,
			'nfit':fit.nfit,
			'bins':len(self.x),
			'dof':len(self.x)-fit.nfit,
			'rchi2':(fit.fval/(len(self.x)-fit.nfit)),
			'fit_xmin':float(fit_xmin),
			'fit_xmax':float(fit_xmax),
			'fit_nsigma':fit_nsigma}
		return par

	def fit_gauss_linear(self,par_init,flag_error=False,fit_nsigma=3,flag_fit_nsigma=True):

		if flag_error:
			chi2reg = Chi2Regression(model_gauss_linear,self.x,self.y,error=self.yerr)
		else:
			chi2reg = Chi2Regression(model_gauss_linear,self.x,self.y)	

		fit = Minuit(chi2reg, 
			peak=par_init['peak'],sigma=par_init['sigma'],area=par_init['area'], 
			c0=par_init['c0'], c1=par_init['c1'],
			limit_peak=(0,None),limit_sigma=(0,None),limit_area=(0,None))
		fit.migrad()
		fit.minos()
		fit.print_param()

		if flag_fit_nsigma:
			peak = fit.values[0]
			sigma = fit.values[1]
			fit_xmin = peak - fit_nsigma * sigma 
			fit_xmax = peak + fit_nsigma * sigma 

			mask = np.logical_and(self.x >= fit_xmin, self.x <= fit_xmax)
			if flag_error:
				chi2reg = Chi2Regression(model_gauss_linear,
					self.x[mask],self.y[mask],error=self.yerr[mask])
			else:
				chi2reg = Chi2Regression(model_gauss_linear,
					self.x[mask],self.y[mask])

			fit = Minuit(chi2reg, 
				peak=par_init['peak'],sigma=par_init['sigma'],area=par_init['area'], 
				c0=par_init['c0'], c1=par_init['c1'])
			fit.migrad()
			fit.minos() 
			fit.print_param()
		else:
			fit_xmin = self.x[0]
			fit_xmax = self.x[-1]

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
			'fval':fit.fval,
			'nfit':fit.nfit,
			'bins':len(self.x),
			'dof':len(self.x)-fit.nfit,
			'rchi2':(fit.fval/(len(self.x)-fit.nfit)),
			'fit_xmin':float(fit_xmin),
			'fit_xmax':float(fit_xmax),
			'fit_nsigma':fit_nsigma}
		return par

	def plot_fit_residual(self,model_x,model_y,outpdf,
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

		mask = (self.yerr > 0.0)

		axs[0] = plt.subplot(gs[0,0])
		if flag_hist:
			axs[0].errorbar(self.x[mask],self.y[mask],
				xerr=self.xerr[mask],yerr=self.yerr[mask],
				marker=None,color='k',drawstyle='steps-mid')
		else:
			axs[0].errorbar(self.x[mask],self.y[mask],
				xerr=self.xerr[mask],yerr=self.yerr[mask],
				marker='o',linestyle='None',color='k')	

		axs[0].plot(model_x,model_y,c='red',drawstyle='steps-mid',linewidth=2)
		axs[0].set_ylabel(ylabel, fontsize=fontsize)
		axs[0].set_xlim(xlim)
		axs[0].set_title(title, fontsize=fontsize)	
		axs[0].get_xaxis().set_visible(False)

		axs[1] = plt.subplot(gs[1])
		if flag_hist:
			axs[1].errorbar(
				self.x[mask],(self.y-model_y)[mask]/self.yerr[mask],
				xerr=self.xerr[mask],yerr=1.0,
				marker=None,linestyle='None',color='k')	
		else:
			axs[1].errorbar(
				self.x[mask],(self.y-model_y)[mask]/self.yerr[mask],
				xerr=self.xerr[mask],yerr=1.0,
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
			if axvline_values is not None:
				for value in axvline_values:
					if axvline_legends is None:
						ax.axvline(value,ls='--')				
					else:
						index = axvline_values.index(value)
						ax.axvline(value,ls='--',label=axvline_legends[index])
				
		legend = axs[0].legend(title=legend_text,loc=legend_loc,fontsize=fontsize_legend)

		#legend.get_title().set_fontsize(fontsize_legend)
		#axs[1].legend(title="(data-model)/error",loc="lower right")

		fig.align_ylabels(axs)
		plt.savefig(outpdf)	

class PhaSpectrum(Hist1D):
	def __init__(self, pha_series, binning, pha_min, pha_max):
		self.pha_min = pha_min
		self.pha_max = pha_max	
		nbins = (self.pha_max - self.pha_min -1)
		nbins = round(nbins / binning)
		super().__init__(nbins, self.pha_min-0.5, self.pha_max-0.5)	
		self.fill(pha_series)

	def write(self,outpdf,title='',xlim=[0,1023]):
		self.plot(outpdf,
			xlabel='ADC channel (pha)',
			ylabel='Counts / bin',
			title=title,
			flag_yerr=True,
			flag_xlog=False,flag_ylog=True,
			xlim=xlim)

class EnergySpectrum(Hist1D):
	def __init__(self, energy_series, nbins, energy_min, energy_max):
		self.energy_min = energy_min
		self.energy_max = energy_max 
		super().__init__(nbins, self.energy_min, self.energy_max)

		self.fill(energy_series)		

	def write(self,outpdf,title='',xlim=[0.0,12.0]):
		self.plot(outpdf,
			xlabel='Energy (MeV)',
			ylabel='Counts / bin',
			title=title,
			flag_yerr=True,
			flag_xlog=False,flag_ylog=True,
			xlim=xlim)

class LightCurve(Hist1D):
	def __init__(self,unixtime_series,unixtime_offset,tbin=8.0,tstart=0.0,tstop=3600.0):
		self.unixtime_series = unixtime_series
		self.unixtime_offset = unixtime_offset 
		self.tbin = tbin
		self.tstart = tstart
		self.tstop = tstop 

		nbins = round((self.tstop - self.tstart)/self.tbin)
		super().__init__(nbins, self.tstart, self.tstop)
		self.fill(self.unixtime_series-self.unixtime_offset)

		self.time_offset_str = datetime.fromtimestamp(self.unixtime_offset)

	def write(self,outpdf,title='',xlim=[0.0,3600.0],mask=None,ylabel=None,
			axvline_values=None,axvline_legends=None,
			axhline_values=None,axhline_legends=None,
			legend_title="",legend_loc="upper right",legend_fontsize=16):
		if ylabel is None:
			ylabel='Counts / (%d sec)' % self.tbin
		self.plot(outpdf,
			xlabel='Time (sec) since %s JST' % self.time_offset_str,
			ylabel=ylabel,
			title=title,
			flag_yerr=True,
			flag_xlog=False,flag_ylog=False,
			xlim=xlim,
			axvline_values=axvline_values,axvline_legends=axvline_values,
			axhline_values=axhline_values,axhline_legends=axhline_values,	
			legend_title=legend_title,legend_loc=legend_loc,legend_fontsize=legend_fontsize,
			mask=mask)

	def set_gaussian_stat(self,outpdf,title='',threshold_sigma=4.0):
		sys.stdout.write('----- {} -----\n'.format(sys._getframe().f_code.co_name))

		ylim = round(1.2*max(self.y))
		self.cnt_hist = Hist1D(nbins=int(ylim),xlow=-0.5,xhigh=ylim-0.5)
		self.cnt_hist.fill(self.y)

		peak = np.mean(self.y)
		sigma = np.std(self.y)
		area = np.sum(self.y)
		fit_xmax = peak + 3 * sigma
		dict_initpar = {'peak':peak,'sigma':sigma,'area':area,'fit_xmin':-0.5,'fit_xmax':fit_xmax}
		print(dict_initpar)
		par = self.cnt_hist.fit_gauss(par_init=dict_initpar,flag_error=None,fit_nsigma=3)		
		par['tbin'] = self.tbin
		par['tstart'] = self.tstart
		par['tstop'] = self.tstop

		model_x = self.cnt_hist.x
		model_y = np.array([model_gauss(x,peak=par['peak'],sigma=par['sigma'],area=par['area']) for x in model_x])	

		self.cnt_hist.plot_fit_residual(model_x,model_y,outpdf,
			xlabel='Count / (%d sec bin)' % self.tbin,
			ylabel='Number of bins',
			title='',flag_hist=True,
			xlim=[-0.5,ylim-0.5])

		return par

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
		self.pdflist = []

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

		sys.stdout.write('[EventData] open {}\n'.format(self.filepath))

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

	def extract_pha_spectrum(self,
			binning=1,pha_min=0,pha_max=MAX_PHA_CHANNEL,
			xlim=[-0.5,MAX_PHA_CHANNEL-0.5]):
		sys.stdout.write('----- {} -----\n'.format(sys._getframe().f_code.co_name))

		self.phaspec = PhaSpectrum(self.df['pha'],
			binning=binning,pha_min=pha_min,pha_max=pha_max)
		outpdf = '%s/%s_phaspec.pdf' % (self.outdir,self.basename)
		title = '%s (%d binning)' % (self.basename,binning)
		self.phaspec.write(outpdf,title=title,xlim=xlim)

		self.pdflist.append(outpdf)
		return outpdf

	def prepare_energy_calibration(self):
		self.line_MeV = []
		self.line_pha = []
		self.line_pha_err = []
		for dict_line in GAMMA_LINES:
			print(dict_line)
			par = self.fit_phaspec_line(
				pha_min=dict_line['pha_min'],
				pha_max=dict_line['pha_max'],
				binning=dict_line['binning'],
				peak=dict_line['peak'],
				sigma=dict_line['sigma'],
				area=dict_line['area'],
				c0=dict_line['c0'],
				c1=dict_line['c1'],
				name=dict_line['name'],
				MeV=dict_line['MeV'],	
				flag_hist=False)	
			self.line_MeV.append(dict_line['MeV'])
			self.line_pha.append(par['peak'])
			self.line_pha_err.append(par['peak_err'])
		self.get_pha2mev_param(
			np.array(self.line_MeV),
			np.array(self.line_pha),np.array(self.line_pha_err))		

	def fit_phaspec_line(self,
			pha_min=202,pha_max=265,binning=2,
			peak=236,sigma=7,area=2600,c0=800,c1=-1,
			flag_hist=False,fit_nsigma=3,name=None,MeV=None):
		sys.stdout.write('----- {} -----\n'.format(sys._getframe().f_code.co_name))

		phaspec = PhaSpectrum(self.df['pha'],
			binning=binning,pha_min=pha_min,pha_max=pha_max)

		peak = phaspec.x[np.argmax(phaspec.y)]
		dict_initpar = {'peak':peak,'sigma':sigma,'area':area,'c0':c0,'c1':c1}		
		par = phaspec.fit_gauss_linear(dict_initpar,flag_error=True,fit_nsigma=fit_nsigma)
		print(par)

		model_x = phaspec.x
		model_y = np.array([model_gauss_linear(x,peak=par['peak'],sigma=par['sigma'],area=par['area'],c0=par['c0'],c1=par['c1']) for x in model_x])	

		if name is None:
			legend_text = 'line fit'
			outpdf = '%s/%s_linefit.pdf' % (self.outdir,self.basename)
		else:
			legend_text = '%s at %.2f MeV (fit range: %d sigma)\n' % (name,MeV,fit_nsigma)
			outpdf = '%s/%s_%s.pdf' % (self.outdir,self.basename,name)
		legend_text += 'peak=%.1f+/-%.1f ch\n' % (par['peak'],par['peak_err'])
		legend_text += 'sigma=%.1f+/-%.1f ch\n' % (par['sigma'],par['sigma_err'])		
		legend_text += 'area=%.1f+/-%.1f counts\n' % (par['area'],par['area_err'])
		legend_text += 'c0=%.1f+/-%.1f counts\n' % (par['c0'],par['c0_err'])		
		legend_text += 'c1=%.1f+/-%.1f counts\n' % (par['c1'],par['c1_err'])			
		legend_text += 'resolution=%.1f %%' % (100.0*get_energy_resolution(par['peak'],par['sigma']))

		fit_range_min = par['peak'] - fit_nsigma * par['sigma']
		fit_range_max = par['peak'] + fit_nsigma * par['sigma']	
		phaspec.plot_fit_residual(
			model_x,model_y,
			outpdf,
			flag_hist=flag_hist,
			xlabel='Channel',ylabel='Counts/bin',title=self.basename,			
			axvline_values=[fit_range_min,fit_range_max],
			axvline_legends=[r"-%d\sigma" % fit_nsigma,r"+%d\sigma" % fit_nsigma],
			legend_text=legend_text,legend_loc='upper right',
			xlim=[par['peak']-1.2*fit_nsigma*par['sigma'],par['peak']+1.2*fit_nsigma*par['sigma']])

		if name is not None:
			self.param[name] = par

		self.pdflist.append(outpdf)	
		return par 	

	def get_pha2mev_param(self,mev_array,pha_array,pha_err_array):
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

		outpdf = '%s/%s_energycal.pdf' % (self.outdir,self.basename)
		plot_xydata(mev_array,pha_array,yerr=pha_err_array,
			model_x=model_x, model_y=model_y,
			outpdf=outpdf,
			xlabel='Energy (MeV)',ylabel='PHA (channel)',title=self.basename)
		self.pdflist.append(outpdf)	

		self.pha2mev_c0 = - mev2pha_c0 / mev2pha_c1
		self.pha2mev_c1 = 1 / mev2pha_c1
		self.param['pha2mev_c0'] = self.pha2mev_c0
		self.param['pha2mev_c1'] = self.pha2mev_c1		
		return self.pha2mev_c0, self.pha2mev_c1

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

	def set_energy_series(self):		
		""" with rand  
		"""
		sys.stdout.write('----- {} -----\n'.format(sys._getframe().f_code.co_name))

		rand_series = np.random.random_sample(len(self.df['pha']))
		rand_pha = self.df['pha'].astype(np.float32) + rand_series - 0.5 
		self.df['energy_mev'] = self.pha2mev_c1 * rand_pha + self.pha2mev_c0

	def extract_energy_spectrum(self,
			nbins=2**10,energy_min=0.0,energy_max=12.0,
			xlim=[0.0,12.0]):
		sys.stdout.write('----- {} -----\n'.format(sys._getframe().f_code.co_name))

		enespec = EnergySpectrum(self.df['energy_mev'],
			nbins=nbins,energy_min=energy_min,energy_max=energy_max)
		outpdf = '%s/%s_energyspec.pdf' % (self.outdir,self.basename)
		title = '%s (%d channel)' % (self.basename,nbins)
		enespec.write(outpdf=outpdf,title=title,xlim=xlim)
		self.pdflist.append(outpdf)			
		return outpdf 

	def extract_curve(self,tbin=8.0,tstart=0.0,tstop=3600.0,
			energy_min=3.0,energy_max=10.0,xlim=[0.0,3600.0]):
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
		lc.write(outpdf=outpdf,title=title,xlim=xlim)	

		lc.mask = mask	
		lc.message = message
		lc.suffix = suffix

		self.pdflist.append(outpdf)			
		return lc

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

	def search_burst(self,lc,threshold_sigma=4.0):
		sys.stdout.write('----- {} -----\n'.format(sys._getframe().f_code.co_name))

		self.par_curve_stat = lc.set_gaussian_stat(
			threshold_sigma=threshold_sigma,
			outpdf='%s/%s_cnthist_gaussfit.pdf' % (self.outdir,self.basename),
			title='%s Burst search (%s) %d-sigma' % (self.basename, lc.message, threshold_sigma)
			)		

		threshold_count = self.par_curve_stat['peak'] + threshold_sigma * self.par_curve_stat['sigma']
		self.par_curve_stat['threshold_sigma'] = threshold_sigma	
		self.par_curve_stat['threshold_count'] = threshold_count

		mask = (lc.y >= threshold_count)
		self.par_curve_stat['burst_time'] = lc.x[mask]

		print(self.par_curve_stat)

		title = '%s (%s) burst candidate %d-sigma' % (self.basename, lc.message, threshold_sigma)
		outpdf = '%s/%s_lc_%s_bstcand.pdf' % (self.outdir,self.basename,lc.suffix)
		lc.write(outpdf=outpdf,title=title,mask=mask)	
		self.pdflist.append(outpdf)					

		gti_start_index = np.argwhere(mask[:-1] < mask[1:])[0] # [False,True] transition
		gti_stop_index = np.argwhere(mask[:-1] > mask[1:])[0] # [False,True] transition

		self.bst_gti_start = np.array([lc.x[gti_start_index]-lc.tbin*0.5])
		self.bst_gti_stop = np.array([lc.x[gti_stop_index]+lc.tbin*0.5])
		self.numof_bst = len(self.bst_gti_start)
		self.bst_list = []

		for i in range(self.numof_bst):
			print(i+1,self.bst_gti_start[i],self.bst_gti_stop[i])
			self.bst_list.append(
				Burst(self,i+1,	self.bst_gti_start[i],self.bst_gti_stop[i]))

	def plot_multi_curves(self,
			tbin=8.0,tstart=0.0,tstop=3600.0,xlim=[0.0,3600.0],
			ebands=[[None,1.0],[1.0,3.0],[3.0,10.0],[10.0,None]]):

		lc_list = []
		for i in range(len(ebands)):
			emin = ebands[i][0]
			emax = ebands[i][1]			
			lc_list.append(self.extract_curve(
				tbin=tbin,tstart=tstart,tstop=tstop,
				energy_min=emin,energy_max=emax,xlim=xlim))

		outpdf = '%s/%s_multiband_lc.pdf' % (self.outdir,self.basename)
		fontsize = 20

		fig, axs = plt.subplots(len(ebands),1, figsize=(8.27,11.69), 
			sharex=True, gridspec_kw={'hspace': 0})
		for i in range(len(ebands)):
			emin = ebands[i][0]
			emax = ebands[i][1]	
			label = '%s-%s MeV' % (emin,emax)
			#axs[i].errorbar(lc_list[i].x,lc_list[i].y,yerr=lc_list[i].yerr,
			#	marker='',drawstyle='steps-mid',color='k')	
			axs[i].step(lc_list[i].x,lc_list[i].y,
				marker='',drawstyle='steps-mid',color='k',label=label)	
			axs[i].set_ylabel(r"Count/(%d sec)" % tbin)
			axs[i].legend(loc="upper right",fontsize=14)
			axs[i].set_xlim(xlim[0],xlim[1])

		axs[-1].set_xlabel('Time (sec) since %s JST' % lc_list[-1].time_offset_str)
		axs[0].set_title(self.basename,fontsize=fontsize)

		for ax in axs:
			ax.label_outer()	
			ax.minorticks_on()
			ax.xaxis.grid(True)
			ax.xaxis.grid(which='major', linestyle='--', color='#000000')
			ax.xaxis.grid(which='minor', linestyle='-.')	
			ax.tick_params(axis="both", which='major', direction='in', length=5)
			ax.tick_params(axis="both", which='minor', direction='in', length=3)	

		fig.align_ylabels(axs)
		plt.tight_layout(pad=2)
		plt.rcParams["font.family"] = "serif"
		plt.rcParams["mathtext.fontset"] = "dejavuserif"	
		plt.savefig(outpdf)

		self.pdflist.append(outpdf)	

		return outpdf

	def analysis_bursts(self,energy_min=3.0,energy_max=10.0,
		tbin=8.0,time_offset=150.0,fit_nsigma=3,tbin_cumlc=1.0):

		mask, message, suffix = self.get_energy_mask(
			energy_min=energy_min,energy_max=energy_max)

		for bst in self.bst_list:
			print("bst-%d" % bst.param["burst_id"])

			tstart=max(bst.param["gti_start"]-time_offset,0)
			tstop=min(bst.param["gti_stop"]+time_offset,3600)
			lc = LightCurve(
				np.array(self.df[mask]['unixtime']),
				float(self.unixtime_offset),
				tbin=tbin,tstart=tstart,tstop=tstop)

			title = '%s (%s) burst-%d' % (self.basename, message, bst.param["burst_id"])
			outpdf = '%s/%s_bst%02d_lc_%s.pdf' % (self.outdir,self.basename,bst.param["burst_id"],suffix)
			lc.write(outpdf=outpdf,title=title,xlim=[tstart,tstop])	
			self.pdflist.append(outpdf)	

			peak = lc.x[np.argmax(lc.y)]
			sigma = (bst.param["gti_stop"] - bst.param["gti_start"])*0.3
			area = sigma * max(lc.y) * 1.3
			c0 = np.mean(lc.y[0:10])

			par = {'peak':peak,'sigma':sigma,'area':area,'c0':c0,'c1':0}
			print(par)

			par = lc.fit_gauss_linear(par,flag_error=True,fit_nsigma=fit_nsigma,flag_fit_nsigma=False)
			print(par)
			par["tbin"] = tbin
			par["tbin_cumlc"] = tbin_cumlc			
			par["energy_min"] = energy_min
			par["energy_max"] = energy_max
			bst.param["lcfit_param"] = par

			model_x = lc.x
			model_y = np.array([model_gauss_linear(x,peak=par['peak'],sigma=par['sigma'],area=par['area'],c0=par['c0'],c1=par['c1']) for x in model_x])	

			outpdf = '%s/%s_bst%02d_lcfit_%s.pdf' % (self.outdir,self.basename,bst.param["burst_id"],suffix)

			legend_text = 'Burst light curve'

			legend_text += 'peak=%.1f+/-%.1f\n' % (par['peak'],par['peak_err'])
			legend_text += 'sigma=%.1f+/-%.1f\n' % (par['sigma'],par['sigma_err'])		
			legend_text += 'area=%.1f+/-%.1f\n' % (par['area'],par['area_err'])
			legend_text += 'c0=%.1f+/-%.1f\n' % (par['c0'],par['c0_err'])		
			legend_text += 'c1=%.1f+/-%.1f\n' % (par['c1'],par['c1_err'])			

			fit_range_min = par['peak'] - fit_nsigma * par['sigma']
			fit_range_max = par['peak'] + fit_nsigma * par['sigma']	
			lc.plot_fit_residual(
				model_x,model_y,
				outpdf,
				flag_hist=True,
				xlabel='Time (sec) since %s JST' % lc.time_offset_str,
				ylabel='Counts / (%d sec)' % lc.tbin,				
				title=self.basename,			
				axvline_values=[fit_range_min,fit_range_max],
				axvline_legends=[r"$-%d\sigma$" % fit_nsigma,r"$+%d\sigma$" % fit_nsigma],
				legend_text=legend_text,legend_loc='upper right',
				xlim=[tstart,tstop]
				)
			self.pdflist.append(outpdf)	

			cumlc = LightCurve(
				np.array(self.df[mask]['unixtime']),
				float(self.unixtime_offset),
				tbin=tbin_cumlc,tstart=tstart,tstop=tstop)
			offset = np.array([model_linear(x,c0=par['c0'],c1=par['c1'])*(tbin_cumlc/tbin) for x in cumlc.x])
			cumlc.y = (cumlc.y.astype(np.float32)-offset).cumsum(dtype=float)

			cumlc.ratio = cumlc.y/cumlc.y[-1]
			at10percent = cumlc.x[cumlc.ratio >= 0.1][0]
			at90percent = cumlc.x[cumlc.ratio >= 0.9][0]
			t80 = at90percent - at10percent
			print(at10percent,at90percent,t80)
			bst.param["t80"] = float(t80)
			bst.param["at10percent"] = float(at10percent)
			bst.param["at90percent"] = float(at90percent)

			title = '%s (%s) burst-%d' % (self.basename, message, bst.param["burst_id"])
			outpdf = '%s/%s_bst%02d_cumlc_%s.pdf' % (self.outdir,self.basename,bst.param["burst_id"],suffix)
			cumlc.write(outpdf=outpdf,title=title,xlim=[tstart,tstop],ylabel='Cumulative (%d sec)' % tbin_cumlc,
				axvline_values=[at10percent,at90percent],
				#axvline_legends=["10 percent","90 percent"],
				axhline_values=[0.1*cumlc.y[-1],0.9*cumlc.y[-1],cumlc.y[-1]],
				#axhline_legends=["10%","90%","100%"],
				legend_loc='upper left'
				)	
			self.pdflist.append(outpdf)	

			bst.set_parameters()
			bst.write_to_yamlfile()

			mask_bst_time = np.logical_and(
				self.df['unixtime'] >= float(self.unixtime_offset) + float(bst.param["at10percent"]),
				self.df['unixtime'] < float(self.unixtime_offset) + float(bst.param["at90percent"]))
			extract_xspec_pha(
				energy_keV_array=np.array(self.df['energy_mev'][mask_bst_time]*1000.0),
				exposure=float(bst.param['t80']),
				outpha='src.pha')

			mask_bgd_time = np.logical_and(
				self.df['unixtime'] >= float(self.unixtime_offset) + float(bst.param["at90percent"] + 100 ),
				self.df['unixtime'] < float(self.unixtime_offset) + float(bst.param["at90percent"] + 500))
			extract_xspec_pha(
				energy_keV_array=np.array(self.df['energy_mev'][mask_bst_time]*1000.0),
				exposure=400.0,
				outpha='bgd.pha')

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

	def pdfmerge(self):
		sys.stdout.write('----- {} -----\n'.format(sys._getframe().f_code.co_name))

		outpdf = '%s/%s_all.pdf' % (self.outdir,self.basename)
		merger = PyPDF2.PdfFileMerger()
		for pdf in self.pdflist:
			merger.append(pdf)
		merger.write(outpdf)
		merger.close()

		return outpdf 

class Burst():
	def __init__(self,eventdata,burst_id,gti_start,gti_stop):
		self.param = {}
		self.eventdata = eventdata
		self.param["burst_id"] = burst_id
		self.param["gti_start"] = float(gti_start)
		self.param["gti_stop"] = float(gti_stop)

		self.param["t80"] = None
		self.param["at10percent"] = None
		self.param["at90percent"] = None
		self.param["lcfit_param"] = None

	def set_parameters(self):	
		self.param["unixtime_peak"] = float(self.eventdata.unixtime_offset) + float(self.param["lcfit_param"]["peak"])
		self.param["jsttime_peak"] = datetime.fromtimestamp(self.param["unixtime_peak"])
		self.param["filename"] = self.eventdata.filename
		self.param["lcfit_param"]["ncount"] = float(self.param["lcfit_param"]["area"])/float(self.param["lcfit_param"]["tbin"])
		self.param["lcfit_param"]["ncount_err"] = float(self.param["lcfit_param"]["area_err"])/float(self.param["lcfit_param"]["tbin"])		

	def write_to_yamlfile(self):
		yamlfile = '%s/%s_bst%02d.yaml' % (self.eventdata.outdir,self.eventdata.basename,self.param["burst_id"])
		with open(yamlfile, "w") as wf:
		    yaml.dump(self.param, wf,default_flow_style=False)		

class Archive(object):
	def __init__(self,parameter_yamlfile):
		self.param = yaml.load(open(parameter_yamlfile),
			Loader=yaml.FullLoader)		
		print("[Archive] %s is generatd." % self.param['archive_name'])

		self.dict = {
			'detid':[],
			'time':[],
			'all':[],
			'allfile':[],
			'lc': [],
			'lcfile': [],
			'pha': [],
			'phafile': [],			
			'spec': [],
			'specfile': [],				
			'process': [],
			'bst': [],
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
			if re.fullmatch(r'\d{3}_\d{8}_\d{2}.csv', os.path.basename(file_path)):
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
		self.dict['all'].append('--')
		self.dict['allfile'].append('--')	
		self.dict['lc'].append('--')
		self.dict['lcfile'].append('--')			
		self.dict['pha'].append('--')
		self.dict['phafile'].append('--')	
		self.dict['spec'].append('--')
		self.dict['specfile'].append('--')						
		self.dict['process'].append('--')
		self.dict['bst'].append('--')
		self.dict['csvlink'].append('<a href="%s">%s</a>' % (file_path,filename))
		self.dict['csvpath'].append(file_path)		
		self.dict['csvfile'].append(filename)
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

		self.df.drop(['csvpath','csvfile','lcfile','phafile','allfile','specfile'],axis=1).to_html('%s/%s.html' % (self.param['outdir'],self.param['archive_name']), render_links=True, escape=False)

	def process(self,index):
		print("[Archive %s] process index of %s" % (self.param['archive_name'],index))

		csvpath = self.df.iloc[index]['csvpath']
		evt = EventData(csvpath)
		outdir = '%s/product/id%s/%s/%s/%s/%s' % (self.param['outdir'],evt.detid_str, evt.year, evt.month, evt.day, evt.hour_jst)
		evt.set_outdir(outdir)

		pdf = evt.extract_pha_spectrum(binning=2)
		self.df.iloc[index]['pha'] = '<a href=\"../%s\">pdf</a>' % (pdf)
		self.df.iloc[index]['phafile'] = pdf

		evt.prepare_energy_calibration()
		evt.set_energy_series()
		evt.set_time_series()
		
		pdf = evt.plot_multi_curves()
		self.df.iloc[index]['lc'] = '<a href=\"../%s\">pdf</a>' % (pdf)
		self.df.iloc[index]['lcfile'] = pdf

		pdf = evt.extract_energy_spectrum()
		self.df.iloc[index]['spec'] = '<a href=\"../%s\">pdf</a>' % (pdf)
		self.df.iloc[index]['specfile'] = pdf

		lc = evt.extract_curve()
		evt.search_burst(lc,threshold_sigma=4.0)
		self.df.iloc[index]['bst'] = '%d' % (len(evt.bst_list))

		if evt.numof_bst > 0:
			evt.analysis_bursts()
		evt.show_summary()
		evt.write_to_fitsfile()
		evt.write_to_yamlfile()	

		pdf = evt.pdfmerge()
		self.df.iloc[index]['all'] = '<a href=\"../%s\">pdf</a>' % (pdf)
		self.df.iloc[index]['allfile'] = pdf

		self.df.iloc[index]['process'] = 'DONE'
