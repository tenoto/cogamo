#!/usr/bin/env python

import os
import argparse

import numpy as np

import cogamo.cogamo as cogamo

__author__ = 'Teruaki Enoto and Miwa Tsurumi'
__version__ = '1.0'
# v1.0, 2022-03-26, the first version

def get_parser():
	"""
	Creates a new argument parser.
	"""
	parser = argparse.ArgumentParser(
		prog="cgm_plot_adc_spectrum.py",
		usage='%(prog)s filelist outcsvfile',
		description="""
extract spectrum from a Cogamo detector event file.
"""	)
	version = '%(prog)s ' + __version__
	parser.add_argument('input_rawcsv', type=str, 
		help='input rawcsv event file.')
	parser.add_argument('-o','--output_pdf', type=str, 
		help='output pdf file.', default=None)
	parser.add_argument('--phabin', type=int, 
		help='input binning (channel).', default=1)	
	parser.add_argument('--phamin', type=int, 
		help='ADC pha minimum to be used for the curve.', default=0)
	parser.add_argument('--phamax', type=int, 
		help='ADC pha maximum to be used for the curve.', default=1023)
	parser.add_argument('--tstart', type=float, 
		help='input time start (tstart) in sec.', default=0.0)
	parser.add_argument('--tstop', type=float, 
		help='input time stop (tstop) in sec.', default=3600.0)	
	parser.add_argument('--flag_rate', type=int, 
		help='0=raw count, 1=rate (divide by time)', default=0)	
	parser.add_argument('--ymin', type=float, 
		help='plot ymin', default=None)
	parser.add_argument('--ymax', type=float, 
		help='plot ymax', default=None)

	return parser

def main(args=None):
	parser = get_parser()
	args = parser.parse_args(args)

	basename = os.path.basename(args.input_rawcsv)
	if args.output_pdf == None:
		args.output_pdf = '%s_phaspec.pdf' % os.path.splitext(basename)[0]

	evtdata = cogamo.EventData(args.input_rawcsv)
	evtdata.set_time_series()
	evtdata.df['offsettime'] = evtdata.df['unixtime'] - float(evtdata.unixtime_offset)

	if args.tstart > 0.0 or args.tstop < 3600.0:
		flag_time = np.logical_and(evtdata.df['offsettime'] >= args.tstart, evtdata.df['offsettime'] <= args.tstop)	
	else:
		flag_time = np.full(len(evtdata.df),True)

	evtdata.phaspec = cogamo.PhaSpectrum(evtdata.df[flag_time]['pha'],
		binning=args.phabin,
		pha_min=args.phamin,
		pha_max=args.phamax)

	title = os.path.basename(basename)
	if (args.ymin is not None) and (args.ymax is not None):
		ylim = [args.ymin,args.ymax]
	else:
		ylim = None 
	
	if args.flag_rate == 1:
		exposure = args.tstop - args.tstart 	
	elif args.flag_rate == 0:
		exposure = 0
	else:
		print("[Error] --flag_rate should be 0 or 1.")
	evtdata.phaspec.write(args.output_pdf,title=title,exposure=exposure,ylim=ylim)


if __name__=="__main__":
	main()