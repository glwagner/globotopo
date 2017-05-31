#!/bin/env python

'''
Access the Smith-Sandwell topography.

Works with versions 8.2 and 9.1b.

'''
import numpy as np
import sys, os, os.path, string

class SSgrid:
    def __init__(self, nlon, nlat, rlt0, rltf):
        self.nlon = nlon
        self.nlat = nlat
        self.rlt0 = rlt0   # Not actually needed
        self.rltf = rltf   # Not actually needed

class SS_file:
    ssgrids = {'8.2': SSgrid(10800, 6336, -72.006, 72.006),
               '9.1b': SSgrid(21600, 17280, -80.738, 80.738)}
    versions = ['9.1b', '8.2'] # highest priority first
    nsubs = [1, 3, 9]  # could go to 21 next

    def __init__(self, nsub=1, version=None):
        '''
        Memory-map a Smith-Sandwell data file.

        nsub is 1, 3, or 9 to select the subsampling;
        version is '8.2' or '9.1b'; if None, the latest
                available will be used.

        The main user method of this object is "extract".
        '''
        sspath = self.find_directory()
        self.sspath = sspath
        if version is None:
            for v in self.versions:
                fname = self.fullpath(sspath, v, nsub)
                if os.path.exists(fname):
                    version = v
                    break
            if version is None:
                raise RuntimeError('Could not find sstopo data file')
        else:
            fname = self.fullpath(sspath, version, nsub)
            if not os.path.exists(fname):
                raise RuntimeError('Could not find file: %s' % fname)
        grid = self.ssgrids[version]
        self.version = version
        self.nsub = nsub
        self.pts_per_degree = (grid.nlon / 360.0) / nsub
        self.nlon = grid.nlon / nsub
        self.nlat = grid.nlat / nsub
        self.fname = fname
        self.z = np.memmap(fname, dtype='>i2', mode='r',
                                shape=(self.nlat, self.nlon))

    def find_directory(self):
        head, tail = os.path.split(os.path.realpath(__file__))
        while 1:
            sspath = os.path.join(head, 'topog', 'sstopo')
            if os.path.isdir(sspath):
                return sspath
            if not tail:
                raise RuntimeError('failed to find topog directory')
            head, tail = os.path.split(head)

    def fullpath(self, sspath, version, nsub):
        if nsub in self.nsubs[1:]:
            fname = os.path.join(sspath, 'topo_%ss%d.img' % (version, nsub))
        elif nsub == 1:
            fname = os.path.join(sspath, 'topo_%s.img' % (version,))
        else:
            raise ValueError('nsub must be in %s' % self.nsubs)
        return fname

    def extract(self, xr, yr, dtype='i2', grid='center'):
        '''
        Return x, y, z for given lon and lat ranges.

        x, y will be 1-D double arrays
        z will be int16 by default; use dtype kwarg for alternatives.

        if grid is 'center' (default), z.shape = (len(y), len(x));
        if grid is 'boundary', then x and y give the cell boundaries,
        and the dimensions of z are reduced by 1.
        '''
        if not grid in ['center', 'boundary']:
            raise RuntimeError('grid must be "center" or "boundary"')
        if grid == 'center':
            ofs = 0.5
        else:
            ofs = 0
        xr = np.asarray(xr)
        yr = np.asarray(yr)
        #If the left side is negative, shift the entire range by 360 to make
        #all longitudes positive.
        if xr[0] < 0:
           xr = xr + 360;
           wrapped = True;
        else:
           wrapped = False;
        ii, jj = self.imgxy2ij(xr, yr[::-1]);
        jj = np.clip(jj, 0, self.nlat-1)
        nlonz = np.diff(ii)[0];   # could check to make sure they are positive
        nlatz = np.diff(jj)[0];

        # Note: z is transposed relative to the input array; Matlab
        # surface and contour functions want Z(lat, lon).
        z = np.zeros((nlatz, nlonz), np.int16);

        # If the right side wraps past 360 degrees, break longitudes into
        # two pieces, the first (ii(1):ii(2)) going up to just befor 360,
        # and the second (ii(3):ii(4)) going from 0 to the right limit.
        if ii[1] > self.nlon:
           ii = np.array([ii[0], self.nlon, 0, ii[1] - self.nlon],
                               dtype=np.int16)

        jjs = np.arange(jj[0], jj[1], dtype=np.int16)
        if grid == 'boundary':
            jjs += 1  # index goes from high lat to low; add 1
                      # to get the index of the lower lat boundary
        ii0 = np.arange(ii[0], ii[1], dtype=np.int16)
        if len(ii) == 2:      # No crossing of the prime meridian; single read.
            x, y = self.imgij2xy(ii0, jjs, ofs=ofs)
            if wrapped: x = x - 360
            z = self.z[slice(*jj), slice(*ii)] # again, here it is (y,x)
        else:                  # Crosses the prime meridian; 2 reads per latitude
            ii1 = np.arange(ii[2], ii[3], dtype=np.int16)
            x, y = self.imgij2xy(np.concatenate((ii0, ii1)), jjs, ofs=ofs);
            nlonz1 = len(ii0)
            nlonz2 = len(ii1)
            z[:, :nlonz1] = self.z[slice(*jj), slice(ii[0], ii[1])]
            z[:, nlonz1:] = self.z[slice(*jj), slice(ii[2], ii[3])]
            if wrapped:  x[:nlonz1] = x[:nlonz1] - 360
                    # x(jz2) was unwrapped when ii was split in two.

        y=np.ascontiguousarray(y[::-1])
        if grid == 'center':
            isl = slice(None, None, -1)
            jsl = slice(None)
        else:
            isl = slice(None, 0, -1)
            jsl = slice(None, nlonz-1)
        z=np.array(z[isl,jsl], dtype=dtype, order='C', copy=True)
        return x, y, z

    def imgxy2ij(self, x, y):
        x = np.asarray(x)
        y = np.asarray(y)
        tanlat = np.tan(y*np.pi/360.0);
        jj = np.array(self.nlat/2 - self.pts_per_degree*180/np.pi *
             np.log( (1+tanlat)/(1-tanlat) ), dtype=np.int16 )
         # -0.5 + self.nlat/2 is the location of the equator in "index space";
         # Note that the grid straddles the equator, because it is symmetric
         # and has an even number of latitudes.
        ii = np.array(x*self.pts_per_degree, dtype=np.int16)
        return ii,jj


    def imgij2xy(self, ii, jj, ofs = 0.5):
        '''
        Calculate lon, lat from indices.

        If ofs is 0.5, these are the centers of the cells;
        if ofs is 0, they are boundaries.
        '''
        ii = np.asarray(ii)
        jj = np.asarray(jj)
        expy = np.exp((jj+ofs-self.nlat/2)*np.pi/(self.pts_per_degree*180))
        y = -2*(180/np.pi)*np.arctan( (expy-1)/(expy+1) )
        x = (ii+ofs)/self.pts_per_degree
        return x, y

    def make_subsets(self, nsubs = None):
        '''
        Generate the subsampled files using a block-median.

        nsubs is a list of odd integer subsamplings; if None,
        it will use the class or instance attribute.

        This will only need to be run when we get a new SS version,
        or if we decide to use other subsamplings.
        '''
        if self.nsub != 1:
            raise RuntimeError("Make subsets only from original file.")
        if nsubs is None:
            nsubs = self.nsubs[1:]
        for nsub in nsubs:
            fname_out = self.fullpath(self.sspath, self.version, nsub)
            fout = open(fname_out, 'wb')
            # temporary array for each new row:
            znew = np.empty((self.nlon/nsub,), dtype='>i2')
            for i0 in range(0, self.nlat, nsub):
                ii = i0/nsub
                for j0 in range(0, self.nlon, nsub):
                    jj = j0/nsub
                    z = self.z[i0:(i0+nsub),j0:(j0+nsub)]
                    znew[jj] = np.median(z, axis=None)
                znew.tofile(fout)
            fout.close()
            # This is very slow; it might be faster to make a
            # normal array from each horizontal strip, and operate
            # on slices of that.

