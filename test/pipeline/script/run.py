#!/usr/bin/env python

import cogamo.cogamo as cogamo

pipe = cogamo.Pipeline('script/parameter.yaml')
pipe.write()

cat = cogamo.Catalog('script/parameter.yaml')
cat.write()

pipe.set_catalog(cat)
#pipe.process(1)
#pipe.write()

for index in range(len(pipe.df)):
	try:
		pipe.process(index)
		pipe.write()
	except:
		continue

exit()

#pipe.write()


"""
#pipe.process(308)
#exit()
for index in range(len(pipe.df)):
	pipe.process(index)
	pipe.write()
"""	