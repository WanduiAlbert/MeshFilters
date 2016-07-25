#!/usr/bin/env python

import numpy as np
import os
import os.path as op
import sys

def loadfirstline(filepath):
    # Read in only the first line
    f = open(filepath, 'r')
    firstline = f.readline()
    f.close()
    return firstline


def traversetree(top):
    ls = []
    walker = os.walk(top)
    for dirpath, dirnames, filenames in walker:
        # print ("I'm entering directory {0}".format(op.basename(dirpath)))
        if filenames: #filenames is not empty
            for afile in filenames:
                afilepath = op.join(dirpath, afile)
                if afile.endswith('txt'):
                    ls += [loadfirstline(afilepath)]
    return ls


if __name__=="__main__":
    if len(sys.argv) == 1:
        print ("Enter the name of the directory to search")
        sys.exit(-1)
    elif len(sys.argv) == 2:
        top = sys.argv[1]
        startlist = traversetree(top)
    else:
        sys.exit(-1)
    startlist = np.array(startlist)
    uniquestarts = np.unique(startlist)
    for entry in uniquestarts:
        print (entry)
    # print (uniquestarts)
    # np.savetxt('uniquefilestarts.txt',uniquestarts)


