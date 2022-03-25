#!/usr/bin/env python

import os
import argparse

import cogamo.cogamo as cogamo

__author__ = 'Teruaki Enoto'
__version__ = '1.0'
# v1.0, 2022-03-26, the first version

def get_parser():
	"""
	Creates a new argument parser.
	"""
	parser = argparse.ArgumentParser(
		prog="cgm_extract_xspecpha.py",
		usage='%(prog)s input_csv',
		description="""
extract xspec pha spectrum from a Cogamo detector rawcsv event file.
"""	)
	version = '%(prog)s ' + __version__
	parser.add_argument('input_rawcsv', type=str, 
		help='input rawcsv event file.')
	parser.add_argument('-o','--output_pha', type=str, 
		help='output pha file name.', default=None)
	return parser

def main(args=None):
	parser = get_parser()
	args = parser.parse_args(args)

	if args.output_pha == None:
		args.output_pha = '%s.pha' % os.path.splitext(os.path.basename(args.input_rawcsv))[0]

	evtdata = cogamo.EventData(args.input_rawcsv)
	cogamo.extract_xspec_pha(
#		evtdata,outpha,exposure,
#	gti=None,start_unixtime=None,stop_unixtime=None):

if __name__=="__main__":
	main()