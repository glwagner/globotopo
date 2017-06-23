import os, sys; sys.path.append('..')
import time
import numpy as np
import matplotlib.pyplot as plt
import globotopo

# Construct a topodata
t1 = time.time()
data = globotopo.SmithSandwell()
print("Elapsed time: {:7.3f}".format(time.time() - t1))

# Test extraction of subsampled global data
glat, glon, gtopo = data.get_all(subsample=64)

# Test extraction of regional data
r1lat, r1lon, r1topo = data.get_region([30,  45, -160, -120], subsample=32)
r2lat, r2lon, r2topo = data.get_region([-20, 40, -30,  20], subsample=32)

# Plot the result
plt.figure()
plt.pcolormesh(glon, glat, gtopo)

plt.figure()
plt.pcolormesh(r1lon, r1lat, r1topo)

plt.figure()
plt.pcolormesh(r2lon, r2lat, r2topo)

plt.show()
