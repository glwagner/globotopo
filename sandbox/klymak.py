# read a smith and sandwell image file: 
# http://topex.ucsd.edu/marine_topo/
from __future__ import division
import numpy as np

def loadtopo():

    nlon = 21600
    nlat = 17280

    fin = open('topo_18.1.img','rb')

    dat = np.fromfile(fin,dtype='i2')
    dat.byteswap('True')

    fin.close()

    # mercator projection!
    rad = np.pi/180.

    arg2 = np.log(np.tan(rad*(45.-80.738/2.)))
    inds = np.arange(nlat)+1.
    arg1 = rad*(nlat-inds+0.5)*1./60.

    lat = 2.*np.arctan(np.exp(arg1+arg2))/rad-90.
    lat = lat[::-1]
    lon = np.linspace(0,360.-1./60.,nlon)+0.5/60.

    dat = dat.reshape([nlat, nlon])
    dat = np.flipud(dat)

    return lat, lon, dat
