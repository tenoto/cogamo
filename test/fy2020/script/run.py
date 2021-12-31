#!/usr/bin/env python

import cogamo.cogamo as cogamo

cogamo.run('script/parameter.yaml')

exit()

pipe = cogamo.Pipeline('script/parameter.yaml')
pipe.write()

cat = cogamo.Catalog('script/parameter.yaml')
cat.write()

pipe.set_catalog(cat)
pipe.process(3)
pipe.write()

exit()

for index in range(len(pipe.df)):
	try:
		pipe.process(index)
		pipe.write()
	except:
		continue


