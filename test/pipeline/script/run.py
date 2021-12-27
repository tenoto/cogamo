#!/usr/bin/env python

import cogamo.cogamo as cogamo

cogamo_archive = cogamo.Archive('script/parameter.yaml')
cogamo_archive.set_csvfiles()
cogamo_archive.convert_to_dataframe()
cogamo_archive.write()
cogamo_archive.process(1)
#cogamo_archive.write()

"""
#cogamo_archive.process(307)
#cogamo_archive.process(308)
#exit()
for index in range(len(cogamo_archive.df)):
	cogamo_archive.process(index)
	cogamo_archive.write()
"""	