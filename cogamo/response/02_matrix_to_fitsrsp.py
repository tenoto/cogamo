#!/usr/bin/env python

import os 
import numpy as np
import pandas as pd
import matplotlib.pylab as plt 
from matplotlib import colors


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
ax.set_ylabel('Output energy (keV)')
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

#### zoom 
img = plt.imshow(H, 
	aspect='equal',
	#cmap=plt.cm.PuRd, 
	cmap=plt.cm.Reds, 
	interpolation='nearest', 
	origin='lower',
	extent=[500,10000,500,10000],
	norm=colors.LogNorm())
fig.savefig('out/02_matrix_zoom_500keV_10MeV.pdf')

img = plt.imshow(H, 
	aspect='equal',
	#cmap=plt.cm.PuRd, 
	cmap=plt.cm.Reds, 
	interpolation='nearest', 
	origin='lower',
	extent=[500,4000,500,4000],
	norm=colors.LogNorm())
fig.savefig('out/02_matrix_zoom_500keV_4MeV.pdf')

img = plt.imshow(H, 
	aspect='equal',
	#cmap=plt.cm.PuRd, 
	cmap=plt.cm.Reds, 
	interpolation='nearest', 
	origin='lower',
	extent=[2000,3000,2000,3000],
	norm=colors.LogNorm())
fig.savefig('out/02_matrix_zoom_2MeV_3MeV.pdf')

