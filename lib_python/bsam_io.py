#!/usr/bin/python
from numpy import *

def parse_adaptive_data_file(filename="out/m00000.dat", maxlevels=100):
    """Returns a data structure that contains all the data from an output file
    written by BSAM."""

    with open(filename, 'r') as f:
        # Read initial data
        time     = float(f.readline())
        maxlevel = int(f.readline())

        # Loop through all levels and patches
        ipatch = 0
        ilevel_current = 0
        levels = []
        this_level = []
        while f.readline():
            ipatch += 1

            # Read initial patch data
            tmp   = f.readline().split()
            ilevel = int(tmp[0])
            ndim  = int(tmp[1])
            r     = int(tmp[2])
            nvars = int(tmp[3])

            # Check if patch is on a new level
            if ilevel != ilevel_current:
                if ilevel > maxlevels:
                    break
                else:
                    levels.append(this_level)
                    this_level = []
                    ilevel_current = ilevel

            # Add new patch to the level array
            this_level.append(parse_patch(f,nvars,ilevel,ipatch))

        # Collect data
        levels.append(this_level)

    return {
        'time'     : time,
        'maxlevel' : maxlevel,
        'npatches' : ipatch,
        'levels'   : levels
    }

def parse_patch(f,nvars,level,ipatch):
    dx   = [float(x) for x in f.readline().split()]
    xmin = [float(x) for x in f.readline().split()]
    xmax = [float(x) for x in f.readline().split()]
    n    = [int(x)   for x in f.readline().split()]
    mg   = [int(x)   for x in f.readline().split()]
    a    = zeros([nvars,n[1]+2, n[0]+2])

    for j in xrange(n[1]+2):
        for i in xrange(n[0]+2):
            a[:,j,i] = [float(x) for x in f.readline().split()]

    return {
        'level'  : level,
        'ipatch' : ipatch,
        'h'      : dx[0],
        'x'      : [xmin[0],xmax[0]],
        'y'      : [xmin[1],xmax[1]],
        'nx'     : n[0],
        'ny'     : n[1],
        'mgx'    : [mg[0],mg[1]],
        'mgy'    : [mg[2],mg[3]],
        'data'   : a
    }

