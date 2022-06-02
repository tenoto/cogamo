#!/usr/bin/env python

import os 
import argparse

import cogamo.cogamo as cogamo

__author__ = 'Teruaki Enoto'
__version__ = '0.02'
# v0.02 : move cogamo_fits to cogamo
# v0.01 : 2020-08-14 : original version

def get_parser():
	"""
	Creates a new argument parser.
	"""
	parser = argparse.ArgumentParser('cgm_plot_hkfile.py',
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description="""
Plot a raw csv-format CoGaMo house keeping data. 
		"""
		)
	version = '%(prog)s ' + __version__
	parser.add_argument('--version', '-v', action='version', version=version,
		help='show version of this command.')
	parser.add_argument('input_hkfile', type=str, 
		help='input rawcsv event file.')
	parser.add_argument('-o','--output_pdf', type=str, 
		help='output pdf file.', default=None)	
	parser.add_argument('--tstart', type=str, 
		help='tstart time (2022-05-30T17:30:00)', default=None)		
	parser.add_argument('--tstop', type=str, 
		help='tstop time (2022-05-30T21:30:00)', default=None)		
	parser.add_argument('--ylog', type=int, 
		help='ylog for the rate (1:log, 0:linear)', default=0)		
	return parser

def main(args=None):
	parser = get_parser()
	args = parser.parse_args(args)
	print(args)

	if args.output_pdf == None:
		args.output_pdf = '%s_hk.pdf' % os.path.splitext(os.path.basename(args.input_hkfile))[0]
	hkdata = cogamo.HKData(args.input_hkfile)
	hkdata.plot(args.output_pdf,tstart=args.tstart,tstop=args.tstop,ylog=args.ylog)

if __name__=="__main__":
	main()