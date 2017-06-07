"""This module defines a 'globotopo' class of topography objects. Key attributes
of this objects are the latitude-longitude grid on which the topography is 
written, and a memory map to the topography data. The globotopo parent class
defines methods that permit the extraction either of all the topography with
globotopo.get_all() or a subset with globotopo.get_region(region)."""

from __future__ import division
import os, time
import numpy as np
from numpy import pi



class topo(object):
    def __init__(self):
        pass

    def get_all(self, subsample=1):
        """Return arrays of latitude, longitude, and topography for the entire
        region spanned by Smith-Sandwell data.
        
        Args:
            subsample (int): Factor by which to subsample the output.
            
        Returns:
            A tuple of numpy arrays that contain latitude, longitude, and 
            topography data, in that order.

        """

        topo = self.topomem[::subsample, ::subsample]
        lat  = self.lat[::subsample]
        lon  = self.lon[::subsample]

        lon, lat = np.meshgrid(lon, lat)

        return lat, lon, topo


    def get_region(self, region, subsample=None):
        """ Extract a rectangular region from the topo data.

        Args:
            region: An array-like input of the form 
                region = [south, north, east, west] that defines the 
                latitude-longitude limits of the region to be extracted.
                The input latitudes must lie between -90 and 90 and the
                input longitudes must lie between -180 and 180.
                For example, to extract a region between 20S and 40N, and
                30W and 5E, set box = [-20, 40, -30, 5].

            subsample (int): Factor by which to subsample the output.

        Returns:
            A tuple of numpy arrays containing the extracted latitude, 
            longitude, and topography data, in that order. 
            
            NOTE: If the region spans the Prime Meridian (lon=0), the western data
            will be assigned negative Longitude values and pasted to the eastern
            data to preserve continuity of the output grid.
        """

        # Extract sides of the box, flipping south/north coordinates if need be
        south, north = np.sort(np.array(region)[[0, 1]])
        west, east = np.array(region)[[2, 3]]

        # Raise hell if something is amiss
        if not -180 <= east <= 180 or not -180 <= west <= 180:
            raise ValueError("Longitudes must lie between +/-180 degrees")
        elif south < self.minlat or north > self.maxlat:
            raise ValueError("Latitudes must lie between +/- "
                    "{} degrees".format(self.maxlat))
        elif east == west or south == north:
            raise ValueError("Latitudes and longitudes must not be unique!")

        # Wrapping is needed if coordinates cross the prime meridian
        if west > east: 
            acrossdateline = True
        else:
            acrossdateline = False

        # Find indices for cutting, taking care not to produce bad indices.
        jsouth = searchsorted_left( self.lat, south)
        jnorth = searchsorted_right(self.lat, north)

        iwest = searchsorted_left( self.lon, west)
        ieast = searchsorted_right(self.lon, east)

        # Cut
        rlat  = self.lat[jsouth:jnorth]

        if not acrossdateline:
            rtopo = self.topomem[jsouth:jnorth, iwest:ieast]
            rlon  = self.lon[iwest:ieast]

        elif acrossdateline:
            nrlat = jnorth - jsouth
            (nreast, nrwest) = (ieast, self.nlon-iwest)

            rlon = np.zeros(nreast+nrwest, dtype=np.float64)
            rtopo = np.zeros((nrlat, nreast+nrwest), dtype=np.float64)

            # Glue together either side of the prime meridian
            rtopo[:, nrwest:] = self.topomem[jsouth:jnorth, :ieast]
            rtopo[:, :nrwest] = self.topomem[jsouth:jnorth, iwest:]

            # Shift coordinates to preserve monotonicity of data
            rlon[nrwest:] = self.lon[:ieast] + 360
            rlon[:nrwest] = self.lon[iwest:]

        if subsample is not None:
            rtopo = rtopo[::subsample, ::subsample]
            rlat  = rlat[::subsample]
            rlon  = rlon[::subsample]

        rlon, rlat = np.meshgrid(rlon, rlat)

        return rlat, rlon, rtopo
       

class smithsandwell(topo):
    def __init__(self, datapath='../data/topo_18.1.img'):
        """Return a smithsandwell globotopo object. The Smith-Sandwell
        grid is a Mercator projection."""

        if not os.path.isfile(datapath):
            raise ValueError( 
                "Check your topo!\n"
                "The file {} does not exist.".format(os.path.abspath(datapath))
            )

        self.datapath = datapath

        # Fixed properties of Smith-Sandwell v18.1
        self.nlon   = 21600
        self.nlat   = 17280
        self.maxlat = 80.738
        self.minlat = -80.738

        # Compute mercator-projected latitude grid
        rad  = pi/180
        arg2 = np.log( np.tan(rad*((90 - self.maxlat)/2)) )
        inds = np.arange(self.nlat)+1.
        arg1 = rad*(self.nlat - inds + 0.5)*1/60

        self.lat = 2*np.arctan( np.exp(arg1 + arg2) )/rad - 90
        self.lat = self.lat[::-1]   # Flips vector over

        # Equaspaced, centered 1/2 arcmin grid in longitude
        self.lon = np.linspace(0, 360-1/60, self.nlon) + 1/120

        # Create memmap of topo data. Coordinates are topomem[lat, lon]
        self.topomem = np.flipud(
            np.memmap(self.datapath, dtype='>i2', shape=(self.nlat, self.nlon), 
            mode='r'))

        # Reshuffte topomem and longitude arrays so longitude lies between +/-180
        imid = self.lon.size//2
        self.lon[imid:] = self.lon[imid:] - 360

        isort = np.argsort(self.lon)

        self.lon = self.lon[isort]
        self.topomem = self.topomem[:, isort]



def searchsorted_left(data, leftside):
    """Return the index of data so that data[ileft] is either the first 
    index in data or lying just left of leftside."""

    if len(data.shape) > 1:
        raise ValueError("Input data must be one-dimensional.")

    ileft = np.max([
        np.searchsorted(data, leftside,  side='left')-1, 0])

    return ileft


def searchsorted_right(data, rightside):
    """Return the index of data so that data[iright] is either the last
    index in data or lying just right of rightside"""

    if len(data.shape) > 1:
        raise ValueError("Input data must be one-dimensional.")

    iright = np.min([
        np.searchsorted(data, rightside,  side='right'), np.size(data)-1])

    return iright


 
