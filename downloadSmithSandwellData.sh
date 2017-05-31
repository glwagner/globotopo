#!/bin/bash

# This script downloads Smith/Sandwell bathymetry data
# from their ftp server.

ftpaddr='ftp://topex.ucsd.edu/pub/global_topo_1min'
version='18.1'
target='./data/'
topofile="topo_$version.img"

# Make target directory if it does not exist
if [ ! -d $target ]; then
    mkdir -p $target
fi

# Use wget to download readme, refs, permissions, and topographic data
wget "$ftpaddr/*.txt" $target
wget "$ftpaddr/$topofile" $target
