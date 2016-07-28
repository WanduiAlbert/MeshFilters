#! /usr/bin/env python

import numpy as np
import os
import os.path as op
import glob
import sys
import re

sample_pattern = re.compile(r'^#( *)SAMPLE\s*$')
nosample_pattern = re.compile(r'^#\s*NO\s*SAMPLE\s*$')
sampleonearm_pattern = re.compile(r'^#( *)SAMPLE( *)?(-+)( *)?ONE')
nosampleonearm_pattern = re.compile(r'^#( *)*NO\s*SAMPLE( *)?(-+)( *)?ONE')
    
def loadfirstline(filepath):
    # Read in only the first line
    f = open(filepath, 'r')
    firstline = f.readline()
    f.close()
    return firstline

def getscans(dirpath, regpattern):
    matchlist = []
    filelist = glob.glob(dirpath + '*.txt')
    for i, filename in enumerate(filelist):
        line1 = loadfirstline(filename)
        if regpattern.search(line1) is not None:
            matchlist += [i]
    return matchlist

def getsampleonearmscans(dirpath):
    return getscans(dirpath, sampleonearm_pattern)

def getnosampleonearmscans(dirpath):
    return getscans(dirpath, nosampleonearm_pattern)

def getsamplescans(dirpath):
    return getscans(dirpath, sample_pattern)

def getnosamplescans(dirpath):
    return getscans(dirpath, nosample_pattern)

if __name__=="__main__":
    if len(sys.argv) == 1:
        print ("Enter the name of the directory to search")
        sys.exit(-1)
    elif len(sys.argv) == 2:
        top = op.join(sys.argv[1], '')
    else:
        sys.exit(-1)

    samplescans = getscans(top, sample_pattern)
    nosamplescans = getscans(top, nosample_pattern)
    sampleonearmscans = getscans(top, sampleonearm_pattern)
    nosampleonearmscans = getscans(top, nosampleonearm_pattern)

    print ("SAMPLE - ONE ARM SCANS: {0}".format(sampleonearmscans))
    print ("NO SAMPLE - ONE ARM SCANS: {0}".format(nosampleonearmscans))
    print ("SAMPLE: {0}".format(samplescans))
    print ("NO SAMPLE: {0}".format(nosamplescans))