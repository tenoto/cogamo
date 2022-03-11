#!/usr/bin/env python

import os 
import sys
import numpy as np
import pandas as pd

from astropy.io import fits

import matplotlib.pylab as plt 
from matplotlib import colors


# Read event list from the Geant4 simulation 
event_list_file = '%s/response_CsI_coated.txt' % os.getenv('COGAMO_RESPONSE_DATA_PATH')
df = pd.read_csv(event_list_file,
	names=["event number","initial energy","deposit energy"],
	delim_whitespace=True)
df['event number'] = df['event number'].astype('int')
total_number = len(df['event number'])
#print(df)
#print(total_number)

"""
edges_input_energy = np.concatenate([
	np.arange(40,400,5),
	np.arange(400,1200,10),
	np.arange(1200,3000,20),
	np.arange(3000,8000,40),
	np.arange(8000,20000,100),
	np.arange(20000,41200,200)])
"""
edges_input_energy = np.arange(40,41200,10000)	
print(edges_input_energy)
print(len(edges_input_energy))

edges_deposit_energy = np.arange(40,41010,20)
print(edges_deposit_energy)
print(len(edges_deposit_energy))

hist2d, xedges, yedges = np.histogram2d(
	df['initial energy'],
	df['deposit energy'],
	bins=(edges_input_energy,edges_deposit_energy))
print(xedges)
print(len(xedges))
print(yedges)
print(len(yedges))

hist1d,bins = np.histogram(df['initial energy'],bins=edges_input_energy)

fig, ax = plt.subplots(1,1, figsize=(11.69,8.27)) # A4 size, inich unit 
fontsize = 18 

bincentres = [(bins[i]+bins[i+1])/2. for i in range(len(bins)-1)]
binwidth = [(-bins[i]+bins[i+1])/2. for i in range(len(bins)-1)]
plt.step(bincentres,hist1d/binwidth,where='mid',color='b',linestyle='--')

#plt.xlabel(xlabel, fontsize=fontsize)
#plt.ylabel(ylabel, fontsize=fontsize)
#plt.title(title,fontsize=fontsize)
#if flag_xlog: plt.xscale('log')				
#if flag_ylog: plt.yscale('log')
#if xlim!=None: plt.xlim(xlim)
	
#ax.minorticks_on()
#ax.grid(True)
#ax.grid(axis='both',which='major', linestyle='--', color='#000000')
#ax.grid(axis='both',which='minor', linestyle='--')	
#ax.tick_params(axis="both", which='major', direction='in', length=5)
#ax.tick_params(axis="both", which='minor', direction='in', length=3)

#plt.tick_params(labelsize=fontsize)
#	plt.rcParams["font.family"] = "serif"
#	plt.rcParams["mathtext.fontset"] = "dejavuserif"	
#	plt.tight_layout()

plt.savefig('hist1d.pdf')

"""
H = hist2d.T
fig = plt.figure(figsize=(8,8),tight_layout=True)
ax = fig.add_subplot(111,title='Cogamo response matrix (original)')
ax.set_xlabel('Input energy (keV)')		
ax.set_ylabel('Deposit energy (keV)')

img = plt.imshow(H, 
	aspect='equal',
	#cmap=plt.cm.PuRd, 
	cmap=plt.cm.Reds, 
	interpolation='nearest', 
	origin='lower',
	extent=[xedges[0],xedges[-1],yedges[0],yedges[-1]],
	norm=colors.LogNorm())
fig.savefig('01_matrix.pdf')
"""

"""
energy_lo_array = np.concatenate([
	np.arange(40,400,5),
	np.arange(400,1200,10),
	np.arange(1200,3000,20),
	np.arange(3000,8000,40),
	np.arange(8000,20000,100),
	np.arange(20000,41000,200)])
energy_hi_array = np.concatenate([energy_lo_array[1:-1],np.array([41000])])

hist2d, xedges, yedges = np.histogram2d(
	df['initial energy'],
	df['deposited energy'],
	bins=(xnbins,ynbins),
	range=((xlow,xhigh),(ylow,yhigh)))
"""