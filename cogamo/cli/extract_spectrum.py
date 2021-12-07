#!/usr/bin/env python

import argparse

import cogamo.cogamo as cogamo

__author__ = 'Teruaki Enoto and Miwa Tsurumi'
__version__ = '0.02'
# v1.0, 2021-12-08, the first version

def get_parser():
	"""
	Creates a new argument parser.
	"""
	parser = argparse.ArgumentParser(
		prog="extract_spectrum.py",
		usage='%(prog)s filelist outcsvfile',
		description="""
extract spectrum from a Cogamo detector event file.
"""	)
	version = '%(prog)s ' + __version__
	parser.add_argument(
		'filelist',metavar='filelist',type=str,
		help='observation file list.') 
#	parser.add_argument(
#		'yamlfile',metavar='yamlfile',type=str,
#		help='parameter file in the yaml format.') 	
	parser.add_argument(
		'outcsvfile',metavar='outcsvfile',type=str,
		help='output csvfile') 		
	return parser

def extract_spectrum(args):
	cogamo.CogamoEvents(args.filelist)
	#csv = xspec.CSVtoXSPEC(args.filelist,args.yamlfile)
	#csv.make_csv2xspec(args.outcsvfile)

def main(args=None):
	parser = get_parser()
	args = parser.parse_args(args)
	extract_spectrum(args)

if __name__=="__main__":
	main()