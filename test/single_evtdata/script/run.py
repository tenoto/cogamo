#!/usr/bin/env python

import cogamo.cogamo as cogamo


pipe = cogamo.Pipeline()
pipe.process_eventdata('data/012_20210108_17.csv')

#evt = cogamo.EventData('data/012_20210108_17.csv')
#evt.extract_pha_spectrum()
#evt.get_enviromental_radiation_lines()
#evt.fit_pha_spectrum_lines()

#evt.generate_process_file()
#evt.show_summary()
#evt.write_to_fitsfile()
