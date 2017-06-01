import os
import sys; sys.path.append('..')
import numpy as np
import matplotlib.pyplot as plt
import globotopo.topodata as topodata

# Construct a topodata
data = topodata.smithsandwell()

# Test extraction of subsampled global data
glat, glon, gtopo = data.get_all(subsample=64)

# Test extraction of regional data
r1lat, r1lon, r1topo = data.get_region([30,  45, 200, 240], subsample=32)
r2lat, r2lon, r2topo = data.get_region([-20, 40, 330,  20], subsample=32)

# Plot the result
plt.figure()
plt.pcolormesh(glon, glat, gtopo)

plt.figure()
plt.pcolormesh(r1lon, r1lat, r1topo)

plt.figure()
plt.pcolormesh(r2lon, r2lat, r2topo)

plt.show()
