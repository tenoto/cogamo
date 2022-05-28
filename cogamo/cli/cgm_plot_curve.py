#!/usr/bin/env python

import os 
import argparse

import numpy as np

import cogamo.cogamo as cogamo

__author__ = 'Teruaki Enoto'
__version__ = '0.01'
# v0.01 : 2022-04-26 : original version

def get_parser():
	"""
	Creates a new argument parser.
	"""
	parser = argparse.ArgumentParser('cgm_plot_curve.py',
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description="""
Plot a raw light curve of the Cogamo detector
		"""
		)
	version = '%(prog)s ' + __version__
	parser.add_argument('--version', '-v', action='version', version=version,
		help='show version of this command.')
	parser.add_argument('input_rawcsv', type=str, 
		help='input rawcsv event file.')
	parser.add_argument('-o','--output_pdf', type=str, 
		help='output pdf file.', default=None)	
	parser.add_argument('--tbin', type=float, 
		help='input time bin (tbin) in sec.', default=8.0)	
	parser.add_argument('--tstart', type=float, 
		help='input time start (tstart) in sec.', default=0.0)
	parser.add_argument('--tstop', type=float, 
		help='input time stop (tstop) in sec.', default=3600.0)	
	parser.add_argument('--phamin', type=int, 
		help='ADC pha minimum to be used for the curve.', default=0)
	parser.add_argument('--phamax', type=int, 
		help='ADC pha maximum to be used for the curve.', default=1023)			
	return parser

def main(args=None):
	parser = get_parser()
	args = parser.parse_args(args)

	basename = os.path.basename(args.input_rawcsv)
	if args.output_pdf == None:
		args.output_pdf = '%s_curve.pdf' % os.path.splitext(basename)[0]

	evtdata = cogamo.EventData(args.input_rawcsv)
	evtdata.set_time_series()
	flag_adc = np.logical_and(evtdata.df['pha'] >= args.phamin, evtdata.df['pha'] <= args.phamax)
	lc = cogamo.LightCurve(np.array(evtdata.df[flag_adc]['unixtime']),
		float(evtdata.unixtime_offset),
		tbin=args.tbin,tstart=args.tstart,tstop=args.tstop)

	title = '%s, %d <= pha <= %d' % (os.path.basename(basename), args.phamin, args.phamax)
	lc.write(outpdf=args.output_pdf,title=title,
		xlim=[args.tstart,args.tstop])	

	#hkdata.plot(args.output_pdf)

if __name__=="__main__":
	main()