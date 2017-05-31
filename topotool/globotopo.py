# read a smith and sandwell image file: 
# http://topex.ucsd.edu/marine_topo/
from __future__ import division
import os
import numpy as np
from numpy import pi

class globaldata(object):

    def __init__(self, topopath='../data/topo_18.1.img'):

        if not os.path.isfile(topopath):
            raise(ValueError, 
                "Check your topo!\n"
                "The file {} does not exist.".format(os.path.abspath(topopath))
            )

        self.topopath = topopath

        self.nlon   = 21600
        self.nlat   = 17280
        self.maxlat = 80.738
        self.minlat = -80.738

        # Mercator projection to get latitudes
        rad  = pi/180
        arg2 = np.log( np.tan(rad*((90 - self.maxlat)/2)) )
        inds = np.arange(self.nlat)+1.
        arg1 = rad*(self.nlat - inds + 0.5)*1/60

        self.lat = 2*np.arctan( np.exp(arg1 + arg2) )/rad - 90
        self.lat = self.lat[::-1]   # Flips vector over

        # Equaspaced, centered 1/2 arcmin grid in longitude
        self.lon = np.linspace(0, 360-1/60, self.nlon) + 1/120

        # Create memmap of topo data
        self.topomem = np.flipud(
            np.memmap(self.topopath, dtype='>i2', shape=(self.nlat, self.nlon))
        )


    def get_global_topo(self, subsample=1):

        with open(self.topopath, 'rb') as topofile:
            topo = np.fromfile(topofile, dtype='i2')
            topo.byteswap('True')

        topo = topo.reshape([self.nlat, self.nlon])
        topo = np.flipud(topo)

        topo = topo[::subsample, ::subsample]
        lat  = self.lat[::subsample]
        lon  = self.lon[::subsample]

        lon, lat = np.meshgrid(lon, lat)

        return lat, lon, topo


    def get_latlon_box(self, box):
        """ Cut a subregion out of memory-mapped topo data.
            'box' is be an array-like input of the form 
            [south, north, east, west]. Note that this code cannot
            fathom with large boxes that cross both the prime meridian 
            and the dateline (yet).
        """

        # Extract sides of the box, flipping south/north coordinates if need be
        south, north = np.sort(np.array(box)[[0, 1]])
        west, east = np.array(box)[[2, 3]]


        # Raise hell if something is amiss
        if east < -180 or west > 180:
            raise(ValueError, 'Longitudes must lie between +/-180 degrees.')
        elif south < self.minlat or north > self.maxlat:
            raise(ValueError, "Latitudes must lie between +/- {} degrees".format(self.maxlat))
        elif east == west or south == north:
            raise(ValueError, "Latitudes and longitudes must not be unique!")

        # Convert to coordinates between 0 and 360
        if west < 0:
            west += 360
        if east < 0:
            east += 360

        # Wrapping is needed if coordinates cross the prime meridian
        if west > east: 
            acrossgreenwich = True
        else:
            acrossgreenwich = False

        # Find indices for cutting
        iwest = np.searchsorted(self.lon, west, side='left')
        ieast = np.searchsorted(self.lon, east, side='right')
        jsouth = np.searchsorted(self.lat, south, side='left')
        jnorth = np.searchsorted(self.lat, north, side='right')

        # Cut
        boxlat  = self.lat[jsouth:jnorth]

        if not acrossgreenwich:

            boxtopo = self.topomem[jsouth:jnorth, iwest:ieast]
            boxlon  = self.lon[iwest:ieast]

        elif acrossgreenwich:

            nboxlat = jnorth - jsouth
            (nboxeast, nboxwest) = (self.nlon-ieast, iwest)
            boxtopo = np.zeros((nboxlat, nboxeast+nboxwest), dtype=np.int16)

            # Glue together either side of the prime meridian
            boxtopo[:, :nboxwest] = self.topomem[jsouth:jnorth, iwest:]
            boxtopo[:, nboxwest:] = self.topomem[jsouth:jnorth, :ieast]

            # Shift coordinates to preserve monotonicity of data
            boxlon[:nboxwest]   = 360-self.lon[ieast:]
            boxlon[nboxeast+1:] = self.lon[:iwest]

        return boxlat, boxlon, boxtopo
