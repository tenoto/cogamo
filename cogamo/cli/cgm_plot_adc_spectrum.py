#!/usr/bin/env python

import os
import argparse

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
	parser.add_argument('--binning', type=int, 
		help='input binning (channel).', default=1)	
	parser.add_argument('--pha_min', type=int, 
		help='input pha_min (channel)', default=0)
	parser.add_argument('--pha_max', type=float, 
		help='input pha_max (channel)', default=cogamo.MAX_PHA_CHANNEL)				
	return parser

def main(args=None):
	parser = get_parser()
	args = parser.parse_args(args)

	if args.output_pdf == None:
		args.output_pdf = '%s_adcspec.pdf' % os.path.splitext(os.path.basename(args.input_rawcsv))[0]
	evtdata = cogamo.EventData(args.input_rawcsv)

	evtdata.phaspec = cogamo.PhaSpectrum(evtdata.df['pha'],
			binning=args.binning,
			pha_min=args.pha_min,
			pha_max=args.pha_max)
	title = os.path.basename(args.input_rawcsv)
	evtdata.phaspec.write(args.output_pdf,title=title)

if __name__=="__main__":
	main()