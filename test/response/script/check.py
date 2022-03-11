#!/usr/bin/env python

from astropy.io import fits

hdu = fits.open('../../cogamo/growth-bgo.rsp')
energy_lo = hdu['MATRIX'].data['ENERG_LO']
energy_hi = hdu['MATRIX'].data['ENERG_HI']
energy_diff = energy_hi - energy_lo
print(energy_diff)
print("5:%d" % sum(energy_diff==5))
print("10:%d" % sum(energy_diff==10))
print("20:%d" % sum(energy_diff==20))
print("40:%d" % sum(energy_diff==40))
print("100:%d" % sum(energy_diff==100))
print("200:%d" % sum(energy_diff==200))

print(hdu['MATRIX'].data['MATRIX'][0])
print(len(hdu['MATRIX'].data['MATRIX'][0]))

