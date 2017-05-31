import os, sys
import numpy as np
import matplotlib.pyplot as plt

# Test if we are in a familiar directory structure or not
upstairs = [ dir for dir in os.walk('../') ]

if ('../topotool' in upstairs and 
    '../data' in upstairs and 
    os.getcwd() is 'tests'):

    sys.path.append('../globotopo/')

    import topotool.globotopo as globotopo
    topopath = '../data/topo_18.1.img'
    
else:
    import globotopo
    topopath = './topo_18.1.img'


# Construct a globotopo object
data = globotopo.globaldata(topopath=topopath)

# Test extraction of subsampled global data
glat, glon, gtopo = data.get_global_topo(subsample=9)

plt.figure()
plt.pcolormesh(glon, glat, gtopo)
plt.show()

# Test extraction of subgrid that does *not* span the Prime Meridian
lat, lon, topo = data.get_latlon_box([1, 25, 120, -160])

plt.figure()
plt.pcolormesh(lon, lat, topo)
plt.show()

# Test extraction of subgrid that spans the Prime Meridian



