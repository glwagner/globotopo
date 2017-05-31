from __future__ import division
import numpy as np
from numpy import pi
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------   
(west, east, south, north) = (10, 40, 30, 40)

topofilepath = ("/data5/glwagner/Software"
    "/topex.ucsd.edu/pub/global_topo_1min/topo_18.1.img")

# Fixed parameters of the Smith-Sandwell dataset
ntopolon    = 21600         # Number of longitude points.
ntopolat    = 17280         # Number of latitude points.
mintopolat  = -80.738       # Most southern extent of grid.
maxtopolat  = 80.738        # Most northern extent of grid.
topodtype   = '>i2'         # Data are big-endian (>) 2 byte signed integers.
nbytes      = 2             # Number of bytes for each datum in file.

arcmin  = 1./60.  # A single arcminute. (1/60 of a degree)
rad     = pi/180. # A single radian.

# Mercator projection transformations. 
# See: https://en.wikipedia.org/wiki/Mercator_projection#Derivation_of_the_Mercator_projection
yR  = lambda phi: np.log( np.tan(pi/4 + phi/2) )
phi = lambda yR:  2*np.arctan(np.exp(yR)) - pi/2

# Construct the Smith and Sandwell grid
topolon = np.arange(0, 360, arcmin)
topolat = phi( np.linspace(yR(maxtopolat*rad), yR(mintopolat*rad), ntopolat) )/rad

# Find the bounding indices of the selected box
ilat = np.searchsorted(topolat, [north, south])
ilon = np.searchsorted(topolon, [west, east])

print(ilat)
print(ilon)

(nlat, nlon) = (ilat[1]-ilat[0], ilon[1]-ilon[0])
(lat, lon) = (topolat[ilat[0]:ilat[1]], topolon[ilon[0]:ilon[1]])

print(lat.shape)
print(lon.shape)

Lat, Lon = np.meshgrid(lat, lon)
bathy = np.ndarray(Lat.shape, dtype='i2')

with open(topofilepath, 'rb') as topofile:
    for i in range(nlat):
        topofile.seek( nbytes*((ilat[0]+i)*ntopolon + ilon[0]) )
        bathy[:, i] = np.fromfile(topofile, dtype=topodtype, count=nlon)

print(Lat.shape)
print(Lon.shape)
print(bathy.shape)

plt.figure()
plt.pcolormesh(Lon, Lat, bathy)
plt.show()
