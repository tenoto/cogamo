#!/usr/bin/env python

import os 
import numpy as np
import pandas as pd
import matplotlib.pylab as plt 
from matplotlib import colors

def crop(image, x1, x2, y1, y2):
	"""
	Return the cropped image at the x1, x2, y1, y2 coordinates
	"""
	if x2 == -1:
		x2=image.shape[1]-1
	if y2 == -1:
		y2=image.shape[0]-1

	mask = np.zeros(image.shape)
	mask[y1:y2+1, x1:x2+1]=1
	m = mask>0

	return image[m].reshape((y2+1-y1, x2+1-x1))

# Read event list from the Geant4 simulation 
event_list_file = '%s/response_CsI_coated.txt' % os.getenv('COGAMO_RESPONSE_DATA_PATH')
df = pd.read_csv(event_list_file,
	names=["event number","initial energy","deposited energy"],
	delim_whitespace=True)
df['event number'] = df['event number'].astype('int')
total_number = len(df['event number'])
print(df)
print(total_number)

xlow   = 40 # keV
xhigh  = 41000 # keV 
xnbins = 2048
#xnbins = 256

ylow   = 40 # keV 
yhigh  = 41000 # keV
ynbins = 2048
#ynbins = xnbins

hist2d, xedges, yedges = np.histogram2d(
	df['initial energy'],
	df['deposited energy'],
	bins=(xnbins,ynbins),
	range=((xlow,xhigh),(ylow,yhigh)))

H = hist2d.T

fig = plt.figure(figsize=(8,8),tight_layout=True)
ax = fig.add_subplot(111,title='Cogamo response matrix (original)')
ax.set_xlabel('Input energy (keV)')		
ax.set_ylabel('Deposited energy (keV)')

img = plt.imshow(H, 
	aspect='equal',
	#cmap=plt.cm.PuRd, 
	cmap=plt.cm.Reds, 
	interpolation='nearest', 
	origin='lower',
	extent=[xedges[0],xedges[-1],yedges[0],yedges[-1]],
	norm=colors.LogNorm())
fig.savefig('01_matrix.pdf')
print(hist2d)

np.savetxt('01_matrix.txt', hist2d, delimiter='\t', fmt='%d')


"""
#### zoom 
H_cropped= crop(H, 0,256,0,256)

img = plt.imshow(H_cropped, 
	aspect='equal',
	#cmap=plt.cm.PuRd, 
	cmap=plt.cm.Reds, 
	#interpolation='nearest', 
	origin='lower',
	#extent=[xedges[0],xedges[-1],yedges[0],yedges[-1]],
	norm=colors.LogNorm())
fig.savefig('01_matrix_zoom.pdf')
"""

