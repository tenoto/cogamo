#!/usr/bin/env python

import os
import argparse

import yaml
import cogamo.cogamo as cogamo

__author__ = 'Teruaki Enoto'
__version__ = '1.0'
# v1.0, 2022-03-26, the first version

def get_parser():
	"""
	Creates a new argument parser.
	"""
	parser = argparse.ArgumentParser(
		prog="cgm_fit_phaspec_line.py",
		usage='%(prog)s input_csv',
		description="""
extract xspec pha spectrum from a Cogamo detector rawcsv event file.
"""	)
	version = '%(prog)s ' + __version__
	parser.add_argument('input_rawcsv', type=str, 
		help='input rawcsv event file.')
	parser.add_argument('--pha_min', type=int, 
		help='input pha_min', default=40)
	parser.add_argument('--pha_max', type=int, 
		help='input pha_max', default=100)	
	parser.add_argument('--binning', type=int, 
		help='input binning', default=2)		
	parser.add_argument('--peak', type=float, 
		help='input peak', default=60)
	parser.add_argument('--sigma', type=float, 
		help='input sigma', default=3)	
	parser.add_argument('--area', type=float, 
		help='input area', default=100)	
	parser.add_argument('--c0', type=float, 
		help='input c0', default=100)	
	parser.add_argument('--c1', type=float, 
		help='input c1', default=-1.0)
	parser.add_argument('--fit_nsigma', type=float, 
		help='input fit_nsigma', default=3)	
	parser.add_argument('--flag_hist', type=bool, 
		help='input flag_hist', default=False)		
	parser.add_argument('--name', type=str, 
		help='input name', default=None)							
	return parser

def main(args=None):
	parser = get_parser()
	args = parser.parse_args(args)

	evtdata = cogamo.EventData(args.input_rawcsv)
	par = evtdata.fit_phaspec_line(
		pha_min=args.pha_min,
		pha_max=args.pha_max,
		binning=args.binning,
		peak=args.peak,
		sigma=args.sigma,
		area=args.area,
		c0=args.c0,
		c1=args.c1,
		flag_hist=args.flag_hist,
		fit_nsigma=args.fit_nsigma,
		name=args.name,
		MeV=None)
	print(par)

	yamlfile = '%s_fit.yaml' % os.path.splitext(os.path.basename(args.input_rawcsv))[0]

	with open(yamlfile, 'w') as outfile:
	    yaml.dump(par, outfile, default_flow_style=True)

if __name__=="__main__":
	main()