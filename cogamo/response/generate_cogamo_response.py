#!/usr/bin/env python

import os 
import sys
import numpy as np
import pandas as pd
import datetime 

from astropy.io import fits

import matplotlib.pylab as plt 
from matplotlib import colors

# http://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/fits/fitsfiles.html
# http://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/docs/memos/cal_gen_92_002/cal_gen_92_002.html#tth_sEc3.2

# ============================
# Set parameter 
VERSION = '0.01'

PARAM_EDGES_INITIAL_ENERGY = [
	[40,400,5],
	[400,1200,10],
	[1200,3000,20],
	[3000,8000,40],
	[8000,20000,100],
	[20000,41200,200]
	] # full 
PARAM_EDGES_DEPOSIT_ENERGY = [40,41000,20] # full 
#PARAM_EDGES_INITIAL_ENERGY = [[40,41000,1000]] # low-resolution 
#PARAM_EDGES_DEPOSIT_ENERGY = [40,41000,500] # low-resolution 
SURFACE_AREA = 75 # cm2 (5x15)
DEGRADATION_SLOPE = -0.50
DEGRADATION_PIVOT = 1.00 # MeV 
DEGRADATION_NORM = 9.22 # percent 
RESP_FILENAME = 'cogamo_fy2020.rsp'
# ============================

product_dir = '%s_prod' % os.path.splitext(RESP_FILENAME)[0]
cmd = 'rm -rf %s; mkdir -p %s' % (product_dir,product_dir) 
os.system(cmd)

# -------------------
# Set Initial energy binning 
tmp_binarray_list = []
for binning in PARAM_EDGES_INITIAL_ENERGY:
	bin_start = binning[0]
	bin_stop = binning[1]
	bin_step = binning[2]
	tmp_binarray_list.append(
		np.arange(bin_start,bin_stop,bin_step))
edges_initial_energy = np.concatenate(tmp_binarray_list)	
nbin_initial_energy = len(edges_initial_energy)-1
#print(edges_initial_energy)

# -------------------
# Set Deposited energ binning 
edges_deposit_energy = np.arange(
	PARAM_EDGES_DEPOSIT_ENERGY[0],
	PARAM_EDGES_DEPOSIT_ENERGY[1]+1,
	PARAM_EDGES_DEPOSIT_ENERGY[2])
nbin_deposit_energy = len(edges_deposit_energy)-1
# print(edges_deposit_energy)

# -------------------
# Read event list from the Geant4 simulation 
# file format: 
# 1: event number
# 2: initial energy of a photon (keV)
# 3: deposit energy to the scintillator (keV)
event_list_file = '%s/response_CsI_coated.txt' % os.getenv('COGAMO_RESPONSE_DATA_PATH')
df = pd.read_csv(event_list_file,
	names=["event number","initial energy","deposit energy"],
	delim_whitespace=True)
df['event number'] = df['event number'].astype('int')
total_number = len(df['event number'])

# -------------------
# Energy resolution 
fwhm = lambda x: x*(DEGRADATION_NORM/100.0)*(x/(DEGRADATION_PIVOT*1000.0))**DEGRADATION_SLOPE if (x>0.0) else 0.0
func_fwhm = np.vectorize(fwhm)
df["fwhm"] = func_fwhm(df["deposit energy"])
df["sigma"] = df["fwhm"] / (2*np.sqrt(2.0*np.log(2.0)))
df["randn"] = np.random.randn(total_number)
df["redist energy"] = df["deposit energy"] + df["randn"] * df["sigma"]

# -------------------
# 1D-Histogram of the initial enegy 
# also get the number of injected photons 
# including escaped non-interacted events 
hist1d_init, _ = np.histogram(
	df['initial energy'],
	bins=edges_initial_energy)
bincenter_initial_energy = (edges_initial_energy[1:]+edges_initial_energy[:-1])/2
binwidth_initial_energy = edges_initial_energy[1:]-edges_initial_energy[:-1]

# -------------------
# [Plot] 1D-Histogram of the initial enegy 
fig, ax = plt.subplots(1,1, figsize=(11.69,8.27)) # A4 size, inich unit 
fontsize = 18 

plt.step(bincenter_initial_energy,hist1d_init/binwidth_initial_energy,
	where='mid',color='r',linestyle='-')
ax.set_xscale('log')
ax.set_xlabel('Initial energy (keV)',fontsize=fontsize)		
ax.set_ylabel('Injected number (photons/bin)',fontsize=fontsize)
ax.set_ylim(0,1.2*max(hist1d_init/binwidth_initial_energy))
ax.tick_params(axis="both", which='major', direction='in', 
	length=5)
ax.tick_params(axis="both", which='minor', direction='in', 
	length=3)
plt.savefig('%s/hist1d_init.pdf' % product_dir)

# -------------------
# 2D-histogram map Count map 
# X: initial energy of a photon (keV)
# Y: deposit energy to the scintillator (keV)
# Z: Number of count 
hist2d_count, _, _ = np.histogram2d(
	df['initial energy'],
	#df['deposit energy'],
	df['redist energy'],
	bins=(edges_initial_energy,edges_deposit_energy))
print("hist2d_count.shape:",hist2d_count.shape)

# -------------------
# normalization 
print("---normalization---")
i = 0
list_normalized_count = []
for i in range(nbin_initial_energy):
	list_normalized_count.append(
		hist2d_count[i]/float(hist1d_init[i]))
hist2d_normalized_count = np.array(list_normalized_count)

# -------------------
# multiply the surface area 
hist2d_resp = SURFACE_AREA * hist2d_normalized_count

# -------------------
# [plot] 2-dim matrix 
fig = plt.figure(figsize=(8,8),tight_layout=True)
ax = fig.add_subplot(111,
	title='Cogamo response matrix (original)')
ax.set_xlabel('Initial energy (keV)')		
ax.set_ylabel('Deposit energy with resolution (keV)')
#ax.set_xscale('log')
#ax.set_yscale('log')
ax.pcolormesh(
	edges_initial_energy,edges_deposit_energy,
	hist2d_count.T,
	shading='flat',rasterized='True',edgecolors='none'
	)
fig.savefig('%s/hist2d_count.pdf' % product_dir)


ax.pcolormesh(
	edges_initial_energy,edges_deposit_energy,
	hist2d_normalized_count.T,
	norm=colors.LogNorm(vmin=0.001,vmax=1.0),
	shading='flat',rasterized='True',edgecolors='none'
	)
fig.savefig('%s/hist2d_normalized_count.pdf' % product_dir)

# -------------
# fits response file parameters
detchans = nbin_deposit_energy + 1 
ngrp = 1
fchan = 1 

# --- prepare EBOUNDS (1st Extension) --- 
channel_array = np.arange(nbin_deposit_energy)
emin_array = edges_deposit_energy[:-1]
emax_array = edges_deposit_energy[1:]
ecen_array = 0.5*(emin_array+emax_array)	

e_min_column = fits.Column(name='E_MIN',format="E",array=emin_array,unit="keV")
e_max_column = fits.Column(name='E_MAX',format="E",array=emax_array,unit="keV")	
channel_column = fits.Column(name='CHANNEL',format='J',array=channel_array,unit="chan")
hdu_ebounds = fits.ColDefs([channel_column,e_min_column,e_max_column])

# --- prepare MATRIX (2nd Extension) --- 
matrix_column = fits.Column(name='MATRIX',
	format="%sE" % nbin_deposit_energy, 
	array=hist2d_resp)
energy_lo_column = fits.Column(name='ENERG_LO',
	format="E",unit="keV",
	array=edges_initial_energy[:-1])
energy_hi_column = fits.Column(name='ENERG_HI',
	format="E",unit="keV",
	array=edges_initial_energy[1:])	
n_grp_column = fits.Column(name='N_GRP',
	format="J",unit="",
	array=np.full(nbin_initial_energy,ngrp))
f_chan_column = fits.Column(name='F_CHAN',
	format="1J",unit="",
	array=np.full(nbin_initial_energy,fchan))
n_chan_column = fits.Column(name='N_CHAN',
	format="1J",unit="",
	array=np.full(nbin_initial_energy,nbin_deposit_energy))
hdu_matrix = fits.ColDefs([energy_lo_column,energy_hi_column,n_grp_column,f_chan_column,n_chan_column,matrix_column])
#hdu_matrix = fits.ColDefs([energy_lo_column,energy_hi_column,n_grp_column,f_chan_column,n_chan_column])

# --- Filling the extensions to a fits file ---
prhdu = fits.PrimaryHDU()
tbhdu_ebounds = fits.BinTableHDU.from_columns(hdu_ebounds)
tbhdu_ebounds.name = 'EBOUNDS'
tbhdu_matrix = fits.BinTableHDU.from_columns(hdu_matrix)
tbhdu_matrix.name = 'MATRIX'	
hdulist = fits.HDUList([prhdu,tbhdu_ebounds,tbhdu_matrix])


tbhdu_ebounds.header['EXTNAME'] = ('EBOUNDS','name of this binary table extension')
tbhdu_ebounds.header['HDUCLAS2'] = ('EBOUNDS','nominal energies of PHA chan boundaries')
tbhdu_ebounds.header['HDUVERS2'] = ('1.3.0   ','Obsolete - included for backwards compatibility')

tbhdu_matrix.header['EXTNAME'] = ('MATRIX','name of this binary table extension')
tbhdu_matrix.header['HDUCLAS2'] = ('RSP_MATRIX','datasets is a spectral response matrix')
tbhdu_matrix.header['HDUVERS2'] = ('1.2.0   ','Obsolete - included for backwards compatibility')
tbhdu_matrix.header['HDUCLAS3'] = ('FULL    ','includes all efficiencies')

for hdu in hdulist:
	hdu.header['HDUCLASS'] = ('OGIP    ','format conforms to OGIP standard')
	hdu.header['HDUCLAS1'] = ('RESPONSE','dataset relates to spectral response')
	hdu.header['HDUVERS']  = ('1.3.0   ','version of format (OGIP memo CAL/GEN/92-002a)')
	hdu.header['HDUVERS1'] = ('1.0.0   ','Obsolete - included for backwards compatibility')
	hdu.header['HDUDOC']  = ('OGIP memos CAL/GEN/92-002 & 92-002a','Documents describing the forma')
	hdu.header['RMFVERSN'] = ('1992a   ','Obsolete - included for backwards compatibility')	

	hdu.header['TELESCOP'] = ('ThundercloudProject','mission/satellite name')
	hdu.header['INSTRUME'] = ('CogamoFY2020','instrument/detector name')
	hdu.header['FILTER'] = ('','filter in use')
	hdu.header['DETNAM'] = ('CsI(Tl):50x50x150mm3','detector in use')
	hdu.header['DETCHANS'] = (detchans,'total number of detector channels')	
	hdu.header['CHANTYPE'] = ('PI','WARNING This is NOT an OGIP-approved value')	
	hdu.header['TLMIN1'] = (channel_array[0],'Minimum value legally allowed in column 1')	
	hdu.header['TLMAX1'] = (channel_array[-1],'Maximum value legally allowed in column 1')	

	hdu.header['COMMENT'] = "-- GROWTH collaboration"
	hdu.header['COMMENT'] = "gennerated by %s (version %s)" % (os.path.basename(__file__),VERSION)
	hdu.header['COMMENT'] = "surface area: 5x15 cm"
	hdu.header['COMMENT'] = "at %s" % datetime.datetime.now()

# --- Write to output fitsfile --- 
cmd = 'rm -f %s' % RESP_FILENAME
print(cmd);os.system(cmd)
hdulist.writeto(RESP_FILENAME)



