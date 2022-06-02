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
import matplotlib.dates as dates
from matplotlib.ticker import MaxNLocator

from iminuit import Minuit
from probfit import Chi2Regression

import PyPDF2

MAX_PHA_CHANNEL = 2**10 - 1 		

FIT_RETRY = 2

#DICT_INITPAR_K40 = {'name':'K40','MeV':1.46083,'peak':132,'sigma':5,'area':18025,'c0':3731,'c1':-21.0,'pha_min':100,'pha_max':184,'binning':2,'xlim':[100,184]}
#DICT_INITPAR_TL208 = {'name':'Tl208','MeV':2.61453,'peak':220,'sigma':7,'area':2651,'c0':798.0,'c1':-3,'pha_min':190,'pha_max':284,'binning':2,'xlim':[190,284],}
#GAMMA_LINES = [DICT_INITPAR_K40,DICT_INITPAR_TL208]

df_Tl208 = pd.read_excel('%s/setenv/energy_calibration_setup.xlsx' % os.getenv('COGAMO_PATH'),sheet_name='Tl208',engine='openpyxl')
df_K40 = pd.read_excel('%s/setenv/energy_calibration_setup.xlsx' % os.getenv('COGAMO_PATH'),sheet_name='K40',engine='openpyxl')
ENERGY_CALIBRATION_SETUP = [df_Tl208,df_K40]

RESPFILE = '%s/cogamo/response/cogamo_fy2020_flat.rsp' % os.getenv('COGAMO_PATH')

plt.rcParams['xtick.major.pad']='8'
plt.rcParams['ytick.major.pad']='8'

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

def extract_xspec_pha(energy_keV_array,outpha,exposure,
	gti=None,start_unixtime=None,stop_unixtime=None):
	number_of_channel = 2048
	energy_min = 40 # keV 
	energy_step = 20 # keV 

	energy_max = energy_min + energy_step * number_of_channel
	hist = Hist1D(nbins=number_of_channel,xlow=energy_min,xhigh=energy_max)
	hist.fill(energy_keV_array)

	f = open('tmp_count.txt','w')
	for i in hist.y:
		f.write("%d\n" % i)
	f.close()

	cmd  = 'ascii2pha infile=tmp_count.txt '
	cmd += 'outfile=%s ' % outpha
	cmd += 'chanpres=no dtype=1 qerror=no rows=- fchan=0 tlmin=0 pois=no filter=NONE '
	cmd += 'detchans=%d ' % number_of_channel
	cmd += 'telescope="ThunderCloud" '
	cmd += 'instrume="CogamoFY2020" '
	cmd += 'detnam="CsI5x15Flat" '
	cmd += 'chantype=PI '
	cmd += 'exposure=%.2f'  % exposure
	print(cmd);os.system(cmd)

	cmd = 'rm -f tmp_count.txt'
	print(cmd);os.system(cmd)

	cmd  = ''	
	if start_unixtime is not None:
		start_time_utc = Time(start_unixtime,format='unix',scale='utc')
		start_time_jst = start_time_utc + timedelta(hours=+9)
		cmd += 'fparkey %.8f %s[1] TSTART comm="Start unixtime in UTC" add=yes;\n' % (start_unixtime,outpha)
		cmd += 'fparkey %s %s[1] STARTJST comm="Start time in JST (ISOT)" add=yes;\n' % (start_time_jst.isot,outpha)
	if stop_unixtime is not None:
		stop_time_utc = Time(stop_unixtime,format='unix',scale='utc')
		stop_time_jst = stop_time_utc + timedelta(hours=+9)			
		cmd += 'fparkey %.8f %s[1] TSTOP comm="Stop unixtime in UTC" add=yes;\n' % (stop_unixtime,outpha)
		cmd += 'fparkey %s %s[1] STOPJST comm="Stop time in JST (ISOT)" add=yes;\n' % (stop_time_jst.isot,outpha)
	print(cmd);os.system(cmd)


def fit_xspec(src_pha,bgd_pha,resp,outdir,basename,
	emin_mev=0.3,emax_mev=40.0,title='',systematic=0.06,
	band1_str='**-0.5,3.0-**',band2_str='**-3.0,10.0-**'):
	sys.stdout.write('----- {} -----\n'.format(sys._getframe().f_code.co_name))

	dump  = "data 1:1 %s\n" % src_pha
	dump += "backgrnd 1 %s\n" % bgd_pha
	dump += "response 1:1 %s\n" % resp
	dump += "mdefine growth_pl E^(-1*index)EXP(-(E/(Ecutoff*1000.0))^alpha) : add\n"
	dump += "model  cflux*growth_pl\n"
	dump += "200      -0.1          0          0      1e+06      1e+06\n"
	dump += "20000    -0.1          0          0      1e+06      1e+06\n"
	dump += "-6       0.01       -100       -100        100        100\n"
	dump += "1.4      0.01        -10		 -10         10         10\n"
	dump += "4.0      0.01       0.01        0.01       100        100\n"
	dump += "1         -1      1e-22      1e-22      1e+22      1e+22\n"
	dump += "1         -1          0          0      1e+20      1e+24\n"

	initial_xcm = "%s/%s_fit_xspec_init.xcm" % (outdir,basename)
	f = open(initial_xcm,'w')
	f.write(dump)
	f.close()

	flog = '%s/%s_fit_xspec_out.log' % (outdir,basename)
	fxcm = '%s/%s_fit_xspec_out.xcm' % (outdir,basename)	
	outps = '%s_fit_xspec_out.ps' % (basename)	

	### plot the spectral file 
	cmd  = 'xspec << EOF\n'
	cmd += '@%s\n' % initial_xcm
	cmd += 'setplot energy mev\n'
	cmd += 'notice **-**\n'
	cmd += 'ignore **-%.1f,%.1f-**\n' % (emin_mev,emax_mev)
	cmd += 'systematic %.4f\n' % systematic
	cmd += 'query yes\n'
	cmd += 'renorm\n'
	cmd += 'fit\n'			
	cmd += 'save all %s\n' % fxcm
	cmd += 'iplot ld del\n'
	cmd += 'csize 1.2\n'
	cmd += 'lab pos y 2.3\n'
	cmd += 'la t %s\n' % title	
	cmd += 'lwid 5\n'
	cmd += 'lwid 5\n'
	cmd += 'lwid 5 on 1..100\n'			
	cmd += 'time off\n'
	cmd += 'lab rotate\n'
	cmd += 'r x 0.2 25\n'
	cmd += 'la y Counts sec\\u-1\\d MeV\\u-1\\d\n'
	cmd += 'col 2 on 2\n'	
	cmd += 'win 2 \n'
	cmd += 'lab 2 col 2 lin 0 100 ls 2 jus lef pos 0.200000003 0 " " \n'
	cmd += 'r y2 -7 7\n'
	cmd += 'lab rotate\n'
	cmd += 'r x 0.2 40.0\n'
	cmd += 'hard %s/cps\n' % outps			
	cmd += 'exit\n'
	cmd += 'log %s\n' % flog
	cmd += 'show all\n'	
	cmd += 'err 1.0 3,4,5\n'
	cmd += 'notice **-**\n'
	cmd += 'ignore %s\n' % band1_str
	cmd += 'show rate\n'
	cmd += 'notice **-**\n'
	cmd += 'ignore %s\n' % band2_str	
	cmd += 'show rate\n'	
	cmd += 'log none\n'		
	cmd += 'exit\n'				
	cmd += 'EOF\n'
	print(cmd)

	fcmd = '%s/%s_fit_xspec.cmd' % (outdir,basename)
	f = open(fcmd,'w')
	f.write(cmd)
	f.close()

	os.system(cmd)

	outpdf = os.path.splitext(outps)[0] + '.pdf'
	os.system('ps2pdf %s' % outps)
	os.system('mv %s %s; rm -f %s' % (outpdf,outdir,outps))

	error_ready = False 
	band1_flag = False
	band2_flag = False	
	for line in open(flog):
		cols = line.split()
		print(error_ready,band1_flag,band2_flag,cols)
		if not len(cols) > 2:
			continue 
		if 'Confidence' in cols:
			error_ready = True
		if cols[1] == '3' and cols[4] == 'lg10Flux':
			lg10Flux = float(cols[6])
		if cols[1] == '4' and cols[4] == 'index':
			index = float(cols[5])	
		if cols[1] == '5' and cols[4] == 'Ecutoff':
			Ecutoff = float(cols[5])
		if error_ready is True and (cols[2] != 'channels'):
			if cols[1] == '3':
				lg10Flux_max = float(cols[3])
				lg10Flux_min = float(cols[2])
				flux = 10**lg10Flux
				flux_err_up = 10**(lg10Flux_max)-10**lg10Flux
				flux_err_down = -10**(lg10Flux_min)+10**lg10Flux
			if cols[1] == '4':
				index_min = float(cols[2])				
				index_max = float(cols[3])
				index_err_up = index_max - index
				index_err_down = index_min - index				
			if cols[1] == '5':
				Ecutoff_min = float(cols[2])	
				Ecutoff_max = float(cols[3])
				Ecutoff_err_up = Ecutoff_min - Ecutoff
				Ecutoff_err_down = Ecutoff_max - Ecutoff				
		if band1_str in str(cols[2]):
			band1_flag = True
		if band1_flag and (str(cols[0]) == '#Net'):
			band1_rate = float(cols[6])
			band1_rate_err = float(cols[8])
			band1_flag = False

		if band2_str in str(cols[2]):
			band2_flag = True
		if band2_flag and (str(cols[0]) == '#Net'):
			band2_rate = float(cols[6])
			band2_rate_err = float(cols[8])
			band2_flag = False

	chisquare = grep(flog,"#Test statistic :",4)[-1]
	dof = grep(flog,"# Null hypothesis probability of",7,dtype=int)[-1]
	reduced_chisquare = chisquare / float(dof)

	fitpar = {}
	fitpar["flux"] = flux
	fitpar["flux_err_up"] = flux_err_up	
	fitpar["flux_err_down"] = flux_err_down		
	fitpar["index"] = index			
	fitpar["index_err_up"] = index_err_up			
	fitpar["index_err_down"] = index_err_down					
	fitpar["Ecutoff"] = Ecutoff				
	fitpar["Ecutoff_err_up"] = Ecutoff_err_up				
	fitpar["Ecutoff_err_down"] = Ecutoff_err_down						
	fitpar["chisquare"] = chisquare	
	fitpar["dof"] = dof	
	fitpar["reduced_chisquare"] = reduced_chisquare	
	fitpar["band1_rate"] = band1_rate	
	fitpar["band1_rate_err"] = band1_rate_err	
	fitpar["band2_rate"] = band2_rate	
	fitpar["band2_rate_err"] = band2_rate_err	
	fitpar["xspec_fit_pdf"]	= '%s/%s' % (outdir,outpdf)
	print(fitpar)

	return fitpar

def grep(logfile,keyword,colnum,dtype=float):
	out_word_list = []
	for line in open(logfile):
		if keyword in line:
			cols = line.split()
			if dtype == float:
				out_word_list.append(float(cols[colnum]))
			elif dtype == int:
				out_word_list.append(int(cols[colnum]))				
	return out_word_list

def run(param_yamlfile):
	pipe = Pipeline(param_yamlfile)
	pipe.write()

	cat = Catalog(param_yamlfile)
	cat.write()

	pipe.set_catalog(cat)
	for index in range(len(pipe.df)):
		try:
			pipe.process(index)
			pipe.write()
		except:
			continue	
	print("DONE!!")

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

	def multiply_factor(self,factor):
		self.y = self.y * factor
		self.yerr = self.yerr * factor

	@property
	def data(self):
		return self.x, self.y

	def plot(self,outpdf,
		xlabel='X title',ylabel='Y title',title='Title',		
		flax_xerr=False,flag_yerr=False,
		flag_xlog=False,flag_ylog=False,
		xlim=None,
		ylim=None,
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
		if ylim is not None: plt.ylim(ylim)		

		ax.minorticks_on()
		ax.grid(True)
		ax.grid(axis='both',which='major', linestyle='--', color='#000000')
		ax.grid(axis='both',which='minor', linestyle='--')	
		ax.tick_params(axis="both", which='major', direction='in', length=5)
		ax.tick_params(axis="both", which='minor', direction='in', length=3)

		if axvline_values is not None:
			if axvline_legends is not None:
				for i in range(len(axvline_values)):
					ax.axvline(axvline_values[i],ls='--',label=axvline_legends[i])
			else:				
				for value in axvline_values:
					ax.axvline(value,ls='--')				

		if axhline_values is not None:
			if axhline_legends is not None:
				for i in range(len(axhline_values)):
					ax.axhline(axhline_values[i],ls='-.',label=axhline_legends[i])
			else:				
				for value in axhline_values:
					ax.axhline(value,ls='-.')				

		if legend_title != "":
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
		for i in range(1,FIT_RETRY+1):
			try:
				fit_valid = fit.migrad()
				if fit_valid:
					break 
			except:
				print("fit error %s" % fit_valid)
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
		for i in range(1,FIT_RETRY+1):
			try:
				fit_valid = fit.migrad()
				if fit_valid:
					break 
			except:
				print("fit error %s" % fit_valid)		
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

	def write(self,outpdf,title='',xlim=[0,1023],ylim=None,exposure=0.0):
		"""
		if exposure == 0, the raw count rate,
		if exposure > 0, the spectrum is divided by the exposure 
		"""
		if exposure == 0.0:
			self.plot(outpdf,
				xlabel='ADC channel (pha)',
				ylabel='Counts / bin',
				title=title,
				flag_yerr=True,
				flag_xlog=False,flag_ylog=True,
				xlim=xlim,ylim=ylim)
		elif exposure > 0.0:
			self.multiply_factor(1.0/exposure)
			self.plot(outpdf,
				xlabel='ADC channel (pha)',
				ylabel='Counts / bin / sec',
				title=title,
				flag_yerr=True,
				flag_xlog=False,flag_ylog=True,
				xlim=xlim,ylim=ylim)			

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

		self.rate = self.y / self.tbin 
		self.rate_err = self.yerr / self.tbin 

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
			axvline_values=axvline_values,axvline_legends=axvline_legends,
			axhline_values=axhline_values,axhline_legends=axhline_legends,	
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
		sys.stdout.write('filepath: %s\n' % self.filepath)

		self.param = {}
		self.pdflist = []

		self.timesys = None		

		self.set_filetype()
		self.open_file()

	def set_filetype(self):
		if re.fullmatch(r'\d{3}_\d{8}_\d{2}.csv', self.filename):
			self.filetype = 'rawcsv'
			self.detid_str, self.yyyymmdd_jst, self.hour_jst = self.basename.split("_")		
			self.year = self.yyyymmdd_jst[0:4]
			self.month = self.yyyymmdd_jst[4:6]
			self.day = self.yyyymmdd_jst[6:8]		
		elif re.match(r'\d{3}_\d{8}_\d{2}', self.filename):
			self.filetype = 'rawcsv'
			self.detid_str, self.yyyymmdd_jst, self.hour_jst, self.note = self.basename.split("_")		
			self.year = self.yyyymmdd_jst[0:4]
			self.month = self.yyyymmdd_jst[4:6]
			self.day = self.yyyymmdd_jst[6:8]			
		else:
			sys.stdout.write("[error] filetype error... at %s\n" % sys._getframe().f_code.co_name)
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
				sys.stdout.write("[error] filetype error... at %s\n" % sys._getframe().f_code.co_name)
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
		for df_line in ENERGY_CALIBRATION_SETUP:
			if int(self.detid_str) in df_line["id"].values:
				line_init = df_line[df_line["id"] == int(self.detid_str)]
			else:
				line_init = df_line[df_line["id"] == 0]		
			par = self.fit_phaspec_line(
				pha_min=float(line_init['pha_min']),
				pha_max=float(line_init['pha_max']),
				binning=int(line_init['binning']),
				peak=float(line_init['peak']),
				sigma=float(line_init['sigma']),
				area=float(line_init['area']),
				c0=float(line_init['c0']),
				c1=float(line_init['c1']),
				name=str(line_init['name'].values[0]),
				MeV=float(line_init['MeV']),	
				flag_hist=False)	
			self.line_MeV.append(float(line_init['MeV']))
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

		if not hasattr(self,'outdir'):
			self.outdir = './'
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
			axvline_legends=[r"-%d$\sigma$" % fit_nsigma,r"+%d$\sigma$" % fit_nsigma],
			legend_text=legend_text,legend_loc='upper right',
			xlim=[par['peak']-1.2*fit_nsigma*par['sigma'],par['peak']+1.2*fit_nsigma*par['sigma']])

		if name is not None:
			self.param[name] = par

		if hasattr(self, 'pdflist'):
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
		self.timesys = 'UTC'

		time_offset_str = str_time + '00:00'
		time_offset_jst = Time(time_offset_str, format='isot', scale='utc', precision=5) 
		self.time_offset_utc = time_offset_jst - timedelta(hours=+9)		
		self.unixtime_offset = self.time_offset_utc.to_value('unix',subfmt='decimal')

	def set_energy_series(self,pha2mev_c1=None,pha2mev_c0=None):		
		""" with rand  
		"""
		sys.stdout.write('----- {} -----\n'.format(sys._getframe().f_code.co_name))

		rand_series = np.random.random_sample(len(self.df['pha']))
		rand_pha = self.df['pha'].astype(np.float32) + rand_series - 0.5 
		if pha2mev_c1 is None:
			pha2mev_c1 = self.pha2mev_c1
		if pha2mev_c0 is None:
			pha2mev_c0 = self.pha2mev_c0
		self.df['energy_mev'] = pha2mev_c1 * rand_pha + pha2mev_c0
		#self.df['energy_mev'] = self.pha2mev_c1 * rand_pha + self.pha2mev_c0

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

	def search_burst(self,lc,threshold_sigma=4.0,catalog=None,maximum_trial=5):
		sys.stdout.write('----- {} -----\n'.format(sys._getframe().f_code.co_name))

		for i in range(maximum_trial):
			try:
				self.par_curve_stat = lc.set_gaussian_stat(
					threshold_sigma=threshold_sigma,
					outpdf='%s/%s_cnthist_gaussfit.pdf' % (self.outdir,self.basename),
					title='%s Burst search (%s) %d-sigma' % (self.basename, lc.message, threshold_sigma))
			except Exception as e:
				sys.stdout.write('Warning %s: trial error %d/%d\n' % (sys._getframe().f_code.co_name, i, maximum_trial))
				pass
			else:
				sys.stdout.write('Log %s: ok %d/%d\n' % (sys._getframe().f_code.co_name, i, maximum_trial))		    	
				break
		else:
			sys.stdout.write('Error %s: trial error %d/%d\n' % (sys._getframe().f_code.co_name, i, maximum_trial))
			pass 

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

		if sum(mask) == 0:
			self.bst_gti_start = None
			self.bst_gti_stop = None
			self.numof_bst = 0
			self.bst_list = None
		else:
			gti_start_index = np.argwhere(mask[:-1] < mask[1:])[0] # [False,True] transition
			gti_stop_index = np.argwhere(mask[:-1] > mask[1:])[0] # [False,True] transition

			self.bst_gti_start = np.array([lc.x[gti_start_index]-lc.tbin*0.5])
			self.bst_gti_stop = np.array([lc.x[gti_stop_index]+lc.tbin*0.5])
			self.numof_bst = len(self.bst_gti_start)
			self.bst_list = []

			for i in range(self.numof_bst):
				print(i+1,self.bst_gti_start[i],self.bst_gti_stop[i])
				self.bst_list.append(
					Burst(self,i+1,	self.bst_gti_start[i],self.bst_gti_stop[i],catalog=catalog))

	def plot_multi_curves(self,
			tbin=8.0,tstart=0.0,tstop=3600.0,xlim=[0.0,3600.0],
			ebands=[[None,1.0],[1.0,3.0],[3.0,10.0],[10.0,None]]):
		sys.stdout.write('----- {} -----\n'.format(sys._getframe().f_code.co_name))	

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

	def extract_xspec_pha(self):
		sys.stdout.write('----- {} -----\n'.format(sys._getframe().f_code.co_name))

		self.xspec_pha = '%s/%s_obs.pha' % (self.outdir,self.basename)
		print('xspec_pha: %s' % self.xspec_pha)
		gti_start_unixtime = float(self.df.head(1)['unixtime'])
		gti_stop_unixtime = float(self.df.tail(1)['unixtime'])
		gti_exposure = gti_stop_unixtime - gti_start_unixtime
		print('GTI: %.6f--%.6f' % (gti_start_unixtime,gti_stop_unixtime))
		print('Exposure: %.3f' % gti_exposure)
		gti = [[gti_start_unixtime],[gti_stop_unixtime]]

		extract_xspec_pha(
			energy_keV_array=np.array(
			self.df['energy_mev'])*1000.0, # convert to keV
			outpha=self.xspec_pha,
			exposure=gti_exposure,
			gti=gti,
			start_unixtime=gti_start_unixtime,
			stop_unixtime=gti_stop_unixtime
			)

	def analysis_bursts(self,energy_min=3.0,energy_max=10.0,
		tbin=8.0,time_offset=150.0,fit_nsigma=3,tbin_cumlc=1.0,
		bgd_tgap_to_src=50,bgd_duration=200):
		sys.stdout.write('----- {} -----\n'.format(sys._getframe().f_code.co_name))

		mask, message, suffix = self.get_energy_mask(
			energy_min=energy_min,energy_max=energy_max)

		for bst in self.bst_list:
			print("bst-%d" % bst.param["burst_id"])

			## generate an energy-filtered lightcurve around the burst candidate
			tstart=max(bst.param["gti_start"]-time_offset,0)
			tstop=min(bst.param["gti_stop"]+time_offset,3600)
			lc = LightCurve(
				np.array(self.df[mask]['unixtime']),
				float(self.unixtime_offset),
				tbin=tbin,tstart=tstart,tstop=tstop)

			## plot the light curve
			title = '%s (%s) burst-%d' % (self.basename, message, bst.param["burst_id"])
			outpdf = '%s/%s_bst%02d_lc_%s.pdf' % (self.outdir,self.basename,bst.param["burst_id"],suffix)
			lc.write(outpdf=outpdf,title=title,xlim=[tstart,tstop])	
			self.pdflist.append(outpdf)	

			## set initial parameters 
			peak = lc.x[np.argmax(lc.y)]
			sigma = (bst.param["gti_stop"] - bst.param["gti_start"])*0.3
			area = sigma * max(lc.y) * 1.3
			c0 = np.mean(lc.y[0:10])

			par = {'peak':peak,'sigma':sigma,'area':area,'c0':c0,'c1':0}
			print("Initial parameters: %s" % par)

			## fit the light curve by a gaussian 
			par = lc.fit_gauss_linear(par,flag_error=True,fit_nsigma=fit_nsigma,flag_fit_nsigma=False)
			print("Gaussian fitting: %s" % par)
			par["tbin"] = tbin
			par["tbin_cumlc"] = tbin_cumlc			
			par["energy_min"] = energy_min
			par["energy_max"] = energy_max
			bst.param["lcfit_param"] = par

			## model Gaussian fitting 
			model_x = lc.x
			model_y = np.array([model_gauss_linear(x,peak=par['peak'],sigma=par['sigma'],area=par['area'],c0=par['c0'],c1=par['c1']) for x in model_x])	

			## plot the fitting light curve
			bst.lcfit_outpdf = '%s/%s_bst%02d_lcfit_%s.pdf' % (self.outdir,self.basename,bst.param["burst_id"],suffix)

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
				bst.lcfit_outpdf,
				flag_hist=True,
				xlabel='Time (sec) since %s JST' % lc.time_offset_str,
				ylabel='Counts / (%d sec)' % lc.tbin,				
				title=self.basename,			
				axvline_values=[fit_range_min,fit_range_max],
				axvline_legends=[r"-%d$\sigma$" % fit_nsigma,r"+%d$\sigma$" % fit_nsigma],
				legend_text=legend_text,legend_loc='upper right',
				xlim=[tstart,tstop]
				)
			self.pdflist.append(bst.lcfit_outpdf)	

			## cummulative light curve
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
			bst.cumlc_outpdf = '%s/%s_bst%02d_cumlc_%s.pdf' % (self.outdir,self.basename,bst.param["burst_id"],suffix)
			cumlc.write(outpdf=bst.cumlc_outpdf,title=title,xlim=[tstart,tstop],ylabel='Cumulative (%d sec)' % tbin_cumlc,
				axvline_values=[at10percent,at90percent],
				axvline_legends=["10 percent at %.1f sec" % at10percent,"90 percent at %.1f sec" % at90percent],
				axhline_values=[0.1*cumlc.y[-1],0.9*cumlc.y[-1],cumlc.y[-1]],
				axhline_legends=["10 percent at %.1f cps" % (0.1*cumlc.y[-1]),"90 percent at %.1f cps" % (0.9*cumlc.y[-1]),"100 percent at %.1f cps" % cumlc.y[-1]],
				legend_loc='lower right'
				)	
			self.pdflist.append(bst.cumlc_outpdf)	

			## extract burst spectral file within T80 
			src_tstart = float(self.unixtime_offset) + float(bst.param["at10percent"])
			src_tstop = float(self.unixtime_offset) + float(bst.param["at90percent"])
			mask_bst_time = np.logical_and(
				self.df['unixtime'] >= src_tstart,
				self.df['unixtime'] < src_tstop)
			src_outpha = '%s/%s_bst%02d_src.pha' % (self.outdir,self.basename,bst.param["burst_id"])
			extract_xspec_pha(
				energy_keV_array=np.array(self.df['energy_mev'][mask_bst_time]*1000.0),
				exposure=float(bst.param['t80']),
				outpha=src_outpha,
				start_unixtime=src_tstart,
				stop_unixtime=src_tstop				
				)

			## extract background spectral file 
			bgd_pre_tstart = max(float(self.unixtime_offset) + float(bst.param["at10percent"]) - bgd_tgap_to_src - bgd_duration,float(self.unixtime_offset))
			bgd_pre_tstop = min(float(self.unixtime_offset) + float(bst.param["at10percent"]) - bgd_tgap_to_src,float(self.unixtime_offset)+3600.0)
			bgd_pre_duration = bgd_pre_tstop - bgd_pre_tstart

			mask_bgd_pre_time = np.logical_and(
				self.df['unixtime'] >= bgd_pre_tstart,
				self.df['unixtime'] < bgd_pre_tstop)
			bgd_pre_outpha = '%s/%s_bst%02d_bgd_pre.pha' % (self.outdir,self.basename,bst.param["burst_id"])
			extract_xspec_pha(
				energy_keV_array=np.array(self.df['energy_mev'][mask_bgd_pre_time]*1000.0),
				exposure=bgd_pre_duration,outpha=bgd_pre_outpha,
				start_unixtime=bgd_pre_tstart,stop_unixtime=bgd_pre_tstop)

			bgd_post_tstart = max(float(self.unixtime_offset) + float(bst.param["at90percent"]) + bgd_tgap_to_src,float(self.unixtime_offset))
			bgd_post_tstop = min(float(self.unixtime_offset) + float(bst.param["at90percent"]) + bgd_tgap_to_src + bgd_duration,float(self.unixtime_offset)+3600.0)
			bgd_post_duration = bgd_post_tstop - bgd_post_tstart

			mask_bgd_post_time = np.logical_and(
				self.df['unixtime'] >= bgd_post_tstart,
				self.df['unixtime'] < bgd_post_tstop)
			bgd_post_outpha = '%s/%s_bst%02d_bgd_post.pha' % (self.outdir,self.basename,bst.param["burst_id"])
			extract_xspec_pha(
				energy_keV_array=np.array(self.df['energy_mev'][mask_bgd_post_time]*1000.0),
				exposure=bgd_post_duration,outpha=bgd_post_outpha,
				start_unixtime=bgd_post_tstart,stop_unixtime=bgd_post_tstop)

			bgd_outpha = '%s/%s_bst%02d_bgd_add.pha' % (self.outdir,self.basename,bst.param["burst_id"])
			cmd = 'ln -s %s .;\n' % bgd_pre_outpha
			cmd += 'ln -s %s .;\n' % bgd_post_outpha			
			cmd += 'mathpha %s+%s ' % (os.path.basename(bgd_pre_outpha),os.path.basename(bgd_post_outpha))
			cmd += 'units=C outfil=%s ' % bgd_outpha
			cmd += 'exposure=CALC areascal=NULL ncomments=0;\n'
			cmd += 'rm -f %s %s;\n' % (os.path.basename(bgd_pre_outpha),os.path.basename(bgd_post_outpha))
			cmd += 'mv %s %s;\n' % (bgd_outpha,self.outdir)
			print(cmd);os.system(cmd)

			### binning spectral file 
			src_bin_outpha = '%s/%s_bst%02d_src_bin.pha' % (self.outdir,self.basename,bst.param["burst_id"])
			script_bin = '%s/%s_bst%02d_grpbin.sh' % (self.outdir,self.basename,bst.param["burst_id"])

			pi_min = 2
			pi_max = 497 # 10 MeV 
			nbin = 26
			f = open(script_bin,'w')
			dump  = '#!/bin/sh -f\n\n'
			dump += 'grppha << EOF\n'
			dump += '%s\n' % src_outpha
			dump += '%s\n' % src_bin_outpha
			tmp_float_pi = np.logspace(
				np.log10(pi_min-1), np.log10(pi_max), nbin, base=10)
			for i in range(len(tmp_float_pi)-1):
				p0 = int(round(tmp_float_pi[i]))+1
				p1 = int(round(tmp_float_pi[i+1]))
				binsize = p1 - p0 + 1 
				dump += 'group %d %d %d\n' % (p0, p1, binsize)
			dump += 'group 498 697 200\n' # 10-14 MeV
			dump += 'group 698 2047 1350\n' # 14-40 MeV
			#dump += 'group 698 997 300\n' # 14-20 MeV
			#dump += 'group 998 1497 500\n' # 20-30 MeV
			#dump += 'group 498 2047 1550\n'
			dump += 'exit\n'
			print(dump)
			f.write(dump)
			f.close()			
			
			cmd = 'chmod +x %s' % script_bin
			os.system(cmd)
			os.system('./%s' % script_bin)

			sub_outps = '%s/%s_bst%02d_sub.ps' % (self.outdir,self.basename,bst.param["burst_id"])
			sub_outpdf = '%s/%s_bst%02d_sub.pdf' % (self.outdir,self.basename,bst.param["burst_id"])

			sub_wbgd_outps = '%s/%s_bst%02d_sub_wbgd.ps' % (self.outdir,self.basename,bst.param["burst_id"])
			sub_wbgd_outpdf = '%s/%s_bst%02d_sub_wbgd.pdf' % (self.outdir,self.basename,bst.param["burst_id"])

			### plot the spectral file 
			cmd  = 'xspec << EOF\n'
			cmd += 'data 1:1 %s\n' % src_bin_outpha
			cmd += 'resp 1:1 %s\n' % RESPFILE
			cmd += 'back 1 %s\n' % bgd_outpha
			cmd += 'setp rebin 0 1\n'
			cmd += 'setplot energy mev\n'
			cmd += 'ignore **-0.2\n'
			cmd += 'iplot ld\n'
			cmd += 'lwid 5\n'
			cmd += 'lwid 5\n'
			cmd += 'lwid 5 on 1..100\n'			
			cmd += 'time off\n'
			cmd += 'la t %s Bst-ID %s\n' % (self.basename,bst.param["burst_id"])
			cmd += 'la y Counts sec\\u-1\\d MeV\\u-1\\d\n'
			cmd += 'lab rotate\n'
			cmd += 'r x 0.2 40.0\n'
			cmd += 'hard %s/cps\n' % sub_outps			
			cmd += 'exit\n'
			cmd += 'exit\n'			
			cmd += 'EOF\n'
			print(cmd)
			fcmd = '%s/%s_bst%02d_sub.xcm' % (self.outdir,self.basename,bst.param["burst_id"])
			f = open(fcmd,'w')
			f.write(cmd)
			f.close()

			os.system(cmd)

			os.system('ps2pdf %s' % sub_outps)
			os.system('mv %s %s' % (os.path.basename(sub_outpdf),self.outdir))
			bst.subspec_pdf = sub_outpdf

			cmd  = 'xspec << EOF\n'
			cmd += 'data 1:1 %s\n' % src_bin_outpha
			cmd += 'resp 1:1 %s\n' % RESPFILE
			cmd += 'back 1 %s\n' % bgd_outpha
			cmd += 'data 2:2 %s\n' % bgd_outpha
			cmd += 'resp 2:2 %s\n' % RESPFILE		
			cmd += 'setp rebin 10 30 1\n'
			cmd += 'setp rebin 0 1 1\n'
			cmd += 'setplot energy mev\n'
			cmd += 'ignore 1-2:**-0.2\n'
			cmd += 'iplot ld\n'
#			cmd += 'line on 2\n'
			cmd += 'lwid 5\n'
			cmd += 'lwid 5\n'
			cmd += 'lwid 5 on 1..100\n'			
			cmd += 'time off\n'
			cmd += 'la t %s Bst-ID %s (w/ bgd)\n' % (self.basename,bst.param["burst_id"])
			cmd += 'la y Counts sec\\u-1\\d MeV\\u-1\\d\n'
			cmd += 'lab rotate\n'
			cmd += 'r x 0.2 40.0\n'
			cmd += 'hard %s/cps\n' % sub_wbgd_outps			
			cmd += 'exit\n'
			cmd += 'exit\n'			
			cmd += 'EOF\n'
			print(cmd)
			fcmd = '%s/%s_bst%02d_sub_wbgd.xcm' % (self.outdir,self.basename,bst.param["burst_id"])
			f = open(fcmd,'w')
			f.write(cmd)
			f.close()
			os.system(cmd)

			os.system('ps2pdf %s' % sub_wbgd_outps)
			os.system('mv %s %s' % (os.path.basename(sub_wbgd_outpdf),self.outdir))
			bst.subspec_wbgd_pdf = sub_wbgd_outpdf

			time_peak_str = datetime.fromtimestamp(float(self.unixtime_offset)+float(par['peak']))

			title  = 'Cogamo ID %s, ' % self.detid_str
			title += '%s (JST), ' % time_peak_str.strftime("%Y-%m-%d %H:%M:%S")
			title += 'T80=%.1f sec (exposure)' % float(bst.param['t80'])
			fit_param = fit_xspec(
				src_pha=src_bin_outpha,
				bgd_pha=bgd_outpha,
				resp=RESPFILE,
				outdir=self.outdir,
				basename='%s_bst%02d' % (self.basename,bst.param["burst_id"]),
				emin_mev=0.3,emax_mev=20.,
				title=title)
			for keyword in fit_param:
				bst.param[keyword] = fit_param[keyword]
			print(bst.param)			

			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			### plot two energy band for publication (summary plot)
			fig, axs = plt.subplots(2,1, figsize=(11.69,8.27), 
				sharex=True, gridspec_kw={'hspace': 0})

			energy_min_threshold = 0.5
			mask1, message1, suffix1 = self.get_energy_mask(
				energy_min=energy_min_threshold,energy_max=energy_min)
			mask2, message2, suffix2 = self.get_energy_mask(
				energy_min=energy_min,energy_max=energy_max)			


			title = 'Cogamo ID %s, ' % self.detid_str
			title += 'time bin %.1f sec, ' % tbin
			title += 'glow %.1f counts (Gauss fit), ' % (par['area']/tbin)
			title += 'T80=%.1f sec' %  bst.param["t80"]
			tstart_summary = -500
			tstop_summary  = 500
			lc1 = LightCurve(
				np.array(self.df[mask1]['unixtime']),
				float(self.unixtime_offset)+par['peak'],
				tbin=tbin,tstart=tstart_summary,tstop=tstop_summary)
			lc2 = LightCurve(
				np.array(self.df[mask2]['unixtime']),
				float(self.unixtime_offset)+par['peak'],
				tbin=tbin,tstart=tstart_summary,tstop=tstop_summary)

			label1 = '%.1f-%.1f MeV'%(energy_min_threshold,energy_min)
			axs[0].errorbar(
				lc1.x,lc1.rate,yerr=lc1.rate_err,fmt='o',label=label1,
				c='k',ls='',markersize=3,linewidth=1,ecolor='k')	
			axs[0].step(
				lc1.x,lc1.rate,
				'-', c='k',mec='k', markersize=2, where='mid')		
			axs[0].set_title(title,fontsize=18)
			axs[0].set_ylabel(r"Count sec$^{-1}$",fontsize=18,labelpad=12)
			axs[0].legend(loc='upper left',fontsize=18)

			#yval_line = 34
			yval_line =0.13* (max(lc1.rate)-min(lc1.rate)) + min(lc1.rate)
			axs[0].plot([bst.param["at10percent"]-par['peak'],bst.param["at90percent"]-par['peak']],[yval_line,yval_line],c='salmon',ls='-',linewidth=5,alpha=0.7)
			axs[0].plot([bst.param["at90percent"]-par['peak']+bgd_tgap_to_src,bst.param["at90percent"]-par['peak']+bgd_tgap_to_src+bgd_duration],[yval_line,yval_line],c='royalblue',ls='-',linewidth=5,alpha=0.7)
			axs[0].plot([bst.param["at10percent"]-par['peak']-bgd_tgap_to_src-bgd_duration,bst.param["at10percent"]-par['peak']-bgd_tgap_to_src],[yval_line,yval_line],c='royalblue',ls='-',linewidth=5,alpha=0.7)			

			label2 = '%.1f-%.1f MeV' % (energy_min,energy_max)
			axs[1].errorbar(
				lc2.x,lc2.rate,yerr=lc2.rate_err,fmt='o',label=label2,
				c='k',ls='',markersize=3,linewidth=1,ecolor='k')	
			axs[1].step(
				lc2.x,lc2.rate,
				'-', c='k',mec='k', markersize=2, where='mid')
			axs[1].set_xlabel("Time (sec) since %s (JST)" % time_peak_str.strftime("%Y-%m-%d %H:%M:%S"), fontsize=18,labelpad=12)
			axs[1].set_ylabel(r"Count sec$^{-1}$",fontsize=18,labelpad=12)			
			axs[1].legend(loc='upper left',fontsize=18)

			for ax in axs:
				ax.minorticks_on()
				#ax.grid(True)
				#ax.grid(axis='both',which='major', linestyle='--', color='#000000')
				#ax.grid(axis='both',which='minor', linestyle='--')	
				ax.tick_params(axis="both",which='major',direction='in',length=6,labelsize=18)
				ax.tick_params(axis="both",which='minor',direction='in',length=6)
				ax.axvline(x=bst.param["at10percent"]-par['peak'],c='r',ls='--')
				ax.axvline(x=bst.param["at90percent"]-par['peak'],c='r',ls='--')
				ax.xaxis.set_major_locator(MaxNLocator(integer=True))

			plt.xticks(fontsize=18) #, rotation=90)
			plt.yticks(fontsize=18) #, rotation=90)
			plt.xlim(tstart_summary,tstop_summary)
			#plt.tick_params(labelsize=18)
			plt.rcParams["font.family"] = "serif"
			plt.rcParams["mathtext.fontset"] = "dejavuserif"	
			plt.tight_layout(pad=2)
			fig.align_ylabels(axs)

			outpdf = '%s/%s_bst%02d_lc.pdf' % (self.outdir,self.basename,bst.param["burst_id"])
			plt.savefig(outpdf)	
			# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

			self.pdflist.append(sub_outpdf)	
			self.pdflist.append(sub_wbgd_outpdf)				
			self.pdflist.append(outpdf)						
			self.pdflist.append(bst.param['xspec_fit_pdf'])				

			bst.set_parameters()
			bst.write_to_yamlfile()
			bst.add_to_catalog()

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

	def run_process(self,outdir,lc_binning=2.0,threshold_sigma=5.0):
		sys.stdout.write('----- {} -----\n'.format(sys._getframe().f_code.co_name))		
		self.set_outdir(outdir)
		self.extract_pha_spectrum(binning=lc_binning)
		self.prepare_energy_calibration()
		self.set_energy_series()
		self.set_time_series()
		self.plot_multi_curves()
		self.extract_energy_spectrum()
		lc = self.extract_curve()
		self.search_burst(lc,threshold_sigma=threshold_sigma)
		if self.numof_bst > 0:
			self.analysis_bursts()
		self.show_summary()
		self.write_to_fitsfile()
		self.write_to_yamlfile()	
		self.pdfmerge()
		print(self.param)

class Burst():
	def __init__(self,eventdata,burst_id,gti_start,gti_stop,catalog=None):
		self.param = {}
		self.eventdata = eventdata
		self.param["burst_id"] = burst_id
		self.param["gti_start"] = float(gti_start)
		self.param["gti_stop"] = float(gti_stop)

		self.param["t80"] = None
		self.param["at10percent"] = None
		self.param["at90percent"] = None
		self.param["lcfit_param"] = None

		self.catalog = catalog 

	def set_catalog(self,catalog):
		self.catalog = catalog

	def set_parameters(self):	
		self.param["unixtime_peak"] = float(self.eventdata.unixtime_offset) + float(self.param["lcfit_param"]["peak"])
		self.param["jsttime_peak"] = datetime.fromtimestamp(self.param["unixtime_peak"])
		self.param["filename"] = self.eventdata.filename
		self.param["lcfit_param"]["ncount"] = float(self.param["lcfit_param"]["area"])/float(self.param["lcfit_param"]["tbin"])
		self.param["lcfit_param"]["ncount_err"] = float(self.param["lcfit_param"]["area_err"])/float(self.param["lcfit_param"]["tbin"])	
		self.param["detid_str"]	= self.eventdata.detid_str

	def write_to_yamlfile(self):
		yamlfile = '%s/%s_bst%02d.yaml' % (self.eventdata.outdir,self.eventdata.basename,self.param["burst_id"])
		with open(yamlfile, "w") as wf:
		    yaml.dump(self.param, wf,default_flow_style=False)		

	def add_to_catalog(self):
		if self.catalog is not None:		
			self.catalog.dict['peaktime'].append(self.param["jsttime_peak"])
			self.catalog.dict['t80'].append(self.param["t80"])
			self.catalog.dict['ncount'].append(self.param["lcfit_param"]["ncount"])	
			self.catalog.dict['ncount_err'].append(self.param["lcfit_param"]["ncount_err"])
			self.catalog.dict['detid'].append(self.param["detid_str"])
			self.catalog.dict['filename'].append(self.param["filename"])
			self.catalog.dict['lc'].append('<a href="../%s">pdf</a>' % self.eventdata.multilc_pdf)	
			self.catalog.dict['spec'].append('<a href="../%s">pdf</a>' % self.subspec_pdf)	
			self.catalog.dict['spec_wbgd'].append('<a href="../%s">pdf</a>' % self.subspec_wbgd_pdf)				
			self.catalog.dict['cumlc'].append('<a href="../%s">pdf</a>' % self.cumlc_outpdf)	
			self.catalog.dict['lcfit'].append('<a href="../%s">pdf</a>' % self.lcfit_outpdf)
			self.catalog.write()

class Pipeline(object):
	def __init__(self,parameter_yamlfile):
		self.param = yaml.load(open(parameter_yamlfile),
			Loader=yaml.FullLoader)		
		print("[Pipeline] %s is generatd." % self.param['pipeline_name'])

		self.catalog = None 

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

		self.set_csvfiles()
		self.convert_dict2df()

	def set_catalog(self,catalog):
		self.catalog = catalog 

	def set_csvfiles(self):
		sys.stdout.write('[Pipeline] {} \n'.format(sys._getframe().f_code.co_name))

		### Find CSV file and make pipelines. 
		csv_filelst = sorted(glob.glob('%s/**/*.csv' % self.param['datadir'],
			recursive=True))
		for file_path in csv_filelst:
			if re.fullmatch(r'\d{3}_\d{8}_\d{2}.csv', os.path.basename(file_path)):
				self.add(file_path)
		print(csv_filelst)	

	def add(self,file_path):
		print("[Pipeline %s] add %s" % (self.param['pipeline_name'],file_path))

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

		if not os.path.exists(self.param['outdir']):
			cmd = 'mkdir -p %s' % self.param['outdir']
			print(cmd);os.system(cmd)

		return 0 

	def convert_dict2df(self):
		self.df = pd.DataFrame.from_dict(self.dict, orient='index').T

	def write(self):

		cmd = 'rm -f %s/%s.{csv,html}' % (self.param['outdir'],self.param['pipeline_name'])
		print(cmd);os.system(cmd)
	

		self.df.to_csv('%s/%s.csv' % (self.param['outdir'],self.param['pipeline_name']))

		self.df.drop(['csvpath','csvfile','lcfile','phafile','allfile','specfile'],axis=1).to_html('%s/%s.html' % (self.param['outdir'],self.param['pipeline_name']), render_links=True, escape=False)

	def process(self,index):
		print("[Pipeline %s] process index of %s" % (self.param['pipeline_name'],index))

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
		
		evt.multilc_pdf = evt.plot_multi_curves()
		self.df.iloc[index]['lc'] = '<a href=\"../%s\">pdf</a>' % (evt.multilc_pdf)
		self.df.iloc[index]['lcfile'] = evt.multilc_pdf

		evt.spec_pdf = evt.extract_energy_spectrum()
		self.df.iloc[index]['spec'] = '<a href=\"../%s\">pdf</a>' % (evt.spec_pdf)
		self.df.iloc[index]['specfile'] = evt.spec_pdf

		lc = evt.extract_curve()
		evt.search_burst(lc,threshold_sigma=self.param['burst_treshold_sigma'],catalog=self.catalog)
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

class Catalog(object):
	def __init__(self,parameter_yamlfile):
		self.param = yaml.load(open(parameter_yamlfile),
			Loader=yaml.FullLoader)		
		print("[Catalog] %s is generatd." % self.param['bstcatalog_name'])

		self.htmlfile = '%s/%s.html' % (self.param['outdir'],self.param['bstcatalog_name'])
		self.csvfile = '%s/%s.csv' % (self.param['outdir'],self.param['bstcatalog_name'])		
		self.outdir = self.param['outdir']

		self.dict = {
			'peaktime':[],	
			't80':[],
			'ncount':[],
			'ncount_err':[],			
			'detid':[],	
			'lc':[],
			'spec':[],
			'spec_wbgd':[],
			'cumlc':[],
			'lcfit':[],
			'filename':[]
			}

		if not os.path.exists(os.path.dirname(self.outdir)):
			cmd = 'mkdir -p %s' % self.outdir
			print(cmd);os.system(cmd)

	def write(self):	
		self.df = pd.DataFrame.from_dict(self.dict, orient='index').T		
		self.df.to_csv(self.csvfile)
		self.df.to_html(self.htmlfile,render_links=True, escape=False)

class HKData():
	"""
	1. yyyy-mm-dd (JST)
	2. HH:MM:SS (JST)
	3. data recording interval (min)
	4. rate1 (cps) below "AREABD1" of " config.csv"
	5. rate2 (cps) between "AREABD1" and  "AREABD2" of " config.csv"
	6. rate3 (cps) between "AREABD2" and  "AREABD3" of " config.csv"
	7. rate4 (cps) between "AREABD3" and  "AREABD4" of " config.csv"
	8. rate5 (cps) between "AREABD4" and  "AREABD5" of " config.csv"
	9. rate6 (cps) above "AREABD5" of " config.csv"
	10. temperature (degC)
	11. pressure (hPa)
	12. humidity (%)
	13. the maximum value among the difference of 10-sec moving average of count rates to the latest count rates (10秒移動平均とCPS値との差の最大値:定義を確認する必要がある) 
	14. optical illumination (lux)
	15. n/a
	16. gps status (0:invalid, 1 or 2: valid)
	17. longitude (deg)
	18. latitude (deg)
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
		self.set_time_series()

	def set_filetype(self):
		if re.fullmatch(r'\d{3}_\d{8}.csv', self.filename):
			self.filetype = 'rawcsv'
			self.detid_str, self.yyyymmdd_jst = self.basename.split("_")		
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
			#cmd = ' nkf -w --overwrite %s' % self.filepath
			#print(cmd);os.system(cmd)

			if self.filetype == 'rawcsv':
				self.df = pd.read_csv(self.filepath, index_col=False, on_bad_lines='skip', encoding='utf-8',
					names=['yyyymmdd','hhmmss','interval','rate1','rate2','rate3','rate4','rate5','rate6','temperature','pressure','humidity','differential','lux','n/a','gps_status','longitude','latitude'],
					dtype={})
				self.nevents = len(self.df)

				# sometimes two lines are collapsed. This process excludes these lines. 
				flag = pd.to_numeric(self.df['latitude'],errors='coerce').isnull()
				self.df = self.df[~flag]

			else:
				sys.stdout.write("[error] filetype error...")
				return -1
		except OSError as e:
			raise

		sys.stdout.write('[HKData] open {}\n'.format(self.filepath))	

	def set_time_series(self):		

		sys.stdout.write('----- {} -----\n'.format(sys._getframe().f_code.co_name))

		tmp_time_series_str = np.char.array(self.df['yyyymmdd'] + 'T' + self.df['hhmmss'])
		tmp_time_series_jst = Time(tmp_time_series_str, format='isot', scale='utc', precision=5) 	
		tmp_time_series_utc = tmp_time_series_jst - timedelta(hours=+9)		
		self.df['unixtime'] = tmp_time_series_utc.to_value('unix',subfmt='decimal')
		self.df['unixtime'] = self.df['unixtime'].astype(np.float64)

	def plot(self,outpdf,tstart=None,tstop=None,ylog=0):
		sys.stdout.write('----- {} -----\n'.format(sys._getframe().f_code.co_name))		
		time_series_utc = Time(self.df['unixtime'],format='unix',scale='utc')
		time_series_jst = time_series_utc.to_datetime(timezone=tz_tokyo)
		plt.rcParams['timezone'] = 'Asia/Tokyo'

		title  = 'DET_ID=%s ' % self.detid_str
		title += '(Longitude=%.3f deg, ' % (np.mean(pd.to_numeric(self.df['latitude'],errors='coerce')))
		title += 'Latitude=%.3f deg)' % (np.mean(pd.to_numeric(self.df['longitude'],errors='coerce')))		
		title += '\n'
		title += '%s ' % str(time_series_jst[0])[0:10]
		title += 'Interval=%d sec ' % (self.df['interval'][0])
		title += '(%s)' % self.filename
		title += '\n'		
		#if self.hdu['HK'].header['AREABD2'] > 0.0:
		#	title += 'Rate L (1+2):<%.1f MeV, ' % (self.hdu['HK'].header['AREABD2']/1000.0)
		#	title += 'Rate M (3+4):%.1f-%.1f MeV, ' % (self.hdu['HK'].header['AREABD2']/1000.0,self.hdu['HK'].header['AREABD4']/1000.0)		
		#	title += 'Rate H (5+6):>%.1f MeV, ' % (self.hdu['HK'].header['AREABD4']/1000.0)
		title += 'Rate L (1+2):<1 MeV, ' 
		title += 'Rate M (3+4):1-3 MeV, '
		title += 'Rate H (5+6):>3 MeV '

		if tstart is not None and tstop is not None:
			tmp_tstart_jst = Time(tstart, format='isot', scale='utc', precision=5) 	
			tmp_tstart_utc = tmp_tstart_jst - timedelta(hours=+9)		
			tstart_jst = tmp_tstart_utc.to_datetime(timezone=tz_tokyo)

			tmp_tstop_jst = Time(tstop, format='isot', scale='utc', precision=5) 	
			tmp_tstop_utc = tmp_tstop_jst - timedelta(hours=+9)		
			tstop_jst = tmp_tstop_utc.to_datetime(timezone=tz_tokyo)		

			flag = np.logical_and(time_series_jst >= tstart_jst,time_series_jst <= tstop_jst)
			time_series_jst = time_series_jst[flag]				
			self.df = self.df[flag]

		# color https://matplotlib.org/stable/tutorials/colors/colors.html
		fig, axs = plt.subplots(8,1, figsize=(8.27,11.69), 
			sharex=True, gridspec_kw={'hspace': 0})
		axs[0].step(
			time_series_jst,
			self.df['rate1']+self.df['rate2'],
			'-', c='salmon',mec='k', markersize=2,where='mid')
		axs[0].set_ylabel(r"Rate L (cps)")
		axs[0].set_title(title)	
		if ylog == 1: axs[0].set_yscale('log')
		axs[1].step(
			time_series_jst,
			self.df['rate3']+self.df['rate4'],
			'-', c='tomato',mec='k', markersize=2,where='mid')
		axs[1].set_ylabel(r"Rate M (cps)")	
		if ylog == 1: axs[1].set_yscale('log')		
		axs[2].step(
			time_series_jst,
			self.df['rate5']+self.df['rate6'],
			'-',c='red',mec='k', markersize=2,where='mid')
		axs[2].set_ylabel(r"Rate H (cps)")				
		if ylog == 1: axs[2].set_yscale('log')				
		axs[3].step(time_series_jst,self.df['temperature'],
			'-',c='k',where='mid')
		axs[3].set_ylabel(r"Temp. (degC)")
		axs[4].step(time_series_jst,self.df['pressure'],
			'-',c='blue',where='mid')
		axs[4].set_ylabel(r"Press. (hPa)")		
		axs[5].step(time_series_jst,self.df['humidity'],
			'-',c='yellowgreen',where='mid')
		axs[5].set_ylabel(r"Humid. (%)")
#		axs[6].step(time_series_jst,self.df['differential'],where='mid')		
#		axs[6].set_ylabel(r"Diff (cps)")	
#		axs[6].set_yscale('log')		
		axs[6].step(time_series_jst,self.df['lux'],
			'-',c='purple',where='mid')		
		axs[6].set_yscale('log')
		axs[6].set_ylabel(r"Illum.")
		axs[7].step(time_series_jst,self.df['gps_status'],
			'-',c='sienna',where='mid')		
		axs[7].set_ylabel(r"GPS status")		
		axs[7].set_xlabel(r"Time (JST)")
		axs[7].set_ylim(-0.5,2.5)
		axs[7].xaxis.set_major_formatter(dates.DateFormatter('%m-%d\n%H:%M'))

		if tstart is not None and tstop is not None:
			tmp_tstart_jst = Time(tstart, format='isot', scale='utc', precision=5) 	
			tmp_tstart_utc = tmp_tstart_jst - timedelta(hours=+9)		
			tstart_jst = tmp_tstart_utc.to_datetime(timezone=tz_tokyo)
			tmp_tstop_jst = Time(tstop, format='isot', scale='utc', precision=5) 	
			tmp_tstop_utc = tmp_tstop_jst - timedelta(hours=+9)		
			tstop_jst = tmp_tstop_utc.to_datetime(timezone=tz_tokyo)
			axs[7].set_xlim(tstart_jst,tstop_jst)
		else:
			axs[7].set_xlim(time_series_jst[0],time_series_jst[-1])

		for ax in axs:
			ax.label_outer()	
			ax.minorticks_on()
			ax.xaxis.grid(True)
			ax.xaxis.grid(which='major', linestyle='--', color='#000000')
			ax.xaxis.grid(which='minor', linestyle='-.')	
			ax.xaxis.set_minor_locator(dates.HourLocator())
			ax.tick_params(axis="both", which='major', direction='in', length=5)
			ax.tick_params(axis="both", which='minor', direction='in', length=3)			

		fig.align_ylabels(axs)
		plt.tight_layout(pad=2)
		plt.rcParams["font.family"] = "serif"
		plt.rcParams["mathtext.fontset"] = "dejavuserif"		
		plt.savefig(outpdf)		
