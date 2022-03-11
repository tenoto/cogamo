#!/usr/bin/env python

import os 
import sys
import numpy as np
import pandas as pd
import matplotlib.pylab as plt 
from matplotlib import colors

from astropy.io import fits

xlow   = 40 # keV
xhigh  = 41000 # keV 
xnbins = 2048

ylow   = 40 # keV 
yhigh  = 41000 # keV
ynbins = 2048

hist2d = np.loadtxt('out/01_matrix.txt')
print(hist2d)

H = hist2d.T
fig = plt.figure(figsize=(8,7),tight_layout=True)
ax = fig.add_subplot(111,title='Cogamo response matrix (reloaded)')
ax.set_xlabel('Input energy (keV)')		
ax.set_ylabel('Deposited energy (keV)')
#if smooth:
#	H=gaussian_filter(H,sigma)
img = plt.imshow(H, 
	aspect='equal',
	#cmap=plt.cm.PuRd, 
	cmap=plt.cm.Reds, 
	interpolation='nearest', 
	origin='lower',
	extent=[xlow,xhigh,ylow,yhigh],
	norm=colors.LogNorm())
fig.savefig('out/02_matrix.pdf')

#NENERGY = 2048
NENERGY = 592
NPI     = 2048
DETCHANS = NPI+1
sys.stdout.write('---- NENERGY=%d, NPI=%d, DETCHANS=%d\n' % (NENERGY,NPI,DETCHANS))

PI_MIN = 40.0
PI_STEP = 20.0
#ENERGY_MIN = 40.0
#ENERGY_STEP = 20.0
N_GRP = 1
F_CHAN = 1

"""
5:72
10:80
20:90
40:125
100:120
200:105
"""

# --- prepare EBOUNDS (1st Extension) --- 
channel_array = np.array([i for i in range(0,NPI)])
e_min_array = np.array([PI_MIN+PI_STEP*i for i in range(0,NPI)])
e_max_array = np.array([PI_MIN+PI_STEP*(i+1) for i in range(0,NPI)])	
e_cen_array = 0.5*(e_min_array+e_max_array)	
e_min_column = fits.Column(name='E_MIN',format="E",array=e_min_array,unit="keV")
e_max_column = fits.Column(name='E_MAX',format="E",array=e_max_array,unit="keV")	
channel_column = fits.Column(name='CHANNEL',format='J',array=channel_array,unit="chan")
hdu_ebounds = fits.ColDefs([channel_column,e_min_column,e_max_column])

##########
rmf_2darray = np.zeros((NENERGY,NPI))

for i_input_energy in hist2d:
	print(len(i_input_energy))
	#for j_deposited_energy in i_input_energy:
		#print(j_deposited_energy)
exit()

energy_lo_array = np.concatenate([
	np.arange(40,400,5),
	np.arange(400,1200,10),
	np.arange(1200,3000,20),
	np.arange(3000,8000,40),
	np.arange(8000,20000,100),
	np.arange(20000,41000,200)])
energy_hi_array = np.concatenate([energy_lo_array[1:-1],np.array([41000])])
print(energy_lo_array)
print(energy_hi_array)

#energy_lo_array = np.array([ENERGY_MIN+ENERGY_STEP*i for i in range(0,NENERGY)])
#energy_hi_array = np.array([ENERGY_MIN+ENERGY_STEP*(i+1) for i in range(0,NENERGY)])	
matrix_form = "%sE" % NPI
matrix_column = fits.Column(name='MATRIX',format=matrix_form, array=rmf_2darray)
#matrix_column = fits.Column(name='MATRIX',format=matrix_form, array=rmf_2darray_norm)
energy_lo_column = fits.Column(name='ENERG_LO',format="E",array=energy_lo_array,unit="keV")
energy_hi_column = fits.Column(name='ENERG_HI',format="E",array=energy_hi_array,unit="keV")	
n_grp_column = fits.Column(name='N_GRP',format="J",array=np.array([N_GRP for i in range(0,NENERGY)]),unit="")
f_chan_column = fits.Column(name='F_CHAN',format="1J",array=np.array([F_CHAN for i in range(0,NENERGY)]),unit="")	
n_chan_column = fits.Column(name='N_CHAN',format="1J",array=np.array([NPI for i in range(0,NENERGY)]),unit="")			
#hdu_matrix = fits.ColDefs([energy_lo_column,energy_hi_column,n_grp_column,f_chan_column,n_chan_column,matrix_column])
hdu_matrix = fits.ColDefs([energy_lo_column,energy_hi_column,n_grp_column,f_chan_column,n_chan_column])


# --- Filling the extensions to a fits file ---
prhdu = fits.PrimaryHDU()
tbhdu_ebounds = fits.BinTableHDU.from_columns(hdu_ebounds)
tbhdu_ebounds.name = 'EBOUNDS'
tbhdu_matrix = fits.BinTableHDU.from_columns(hdu_matrix)
tbhdu_matrix.name = 'MATRIX'	
hdu = fits.HDUList([prhdu,tbhdu_ebounds,tbhdu_matrix])
#hdu = fits.HDUList([prhdu,tbhdu_ebounds])


# --- Write to output fitsfile --- 
cmd = 'rm -f test.rsp'
print(cmd);os.system(cmd)
hdu.writeto('test.rsp')


