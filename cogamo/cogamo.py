# -*- coding: utf-8 -*-

class CogamoEvents():
	"""
	Event list of cogamo detectors imported from a csv file (files).

	Input csvfile: 
		-- The event data file recorded all the individual radiation events.  
		-- The file name is [det_id]_[yyyymmdd]_[hour].csv (e.g., 012_20210108_17.csv)
		-- The time is defined in JST.
	The event data file include the following columns:
		1. minute
		2. sec
		3. 1/10000 sec
		4. ADC channel (pha: pulse height amplitude) [0-1023]
	Example: 
		0,0,77,60
		0,0,122,23
		0,0,166,41
	"""
	def __init__(self,csvfilename):
		print("csvfilename=%s" % csvfilename)

	def writeToCSVfile(self,outcsvfile):
		print("outcsvfile=%s" % outcsvfile)

class CogamoSpectrum():
	"""
	Spectrum
	"""
	def __init__(self,csvfilename):
		print("csvfilename=%s" % csvfilename)

class CogamoCurve():
	"""
	Light curve
	"""
	def __init__(self,csvfilename):
		print("csvfilename=%s" % csvfilename)

