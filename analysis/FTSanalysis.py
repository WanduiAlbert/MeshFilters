#!/usr/bin/env python

import numpy as np
import FTSscan as fs
import sys
import os.path as op
import os
import sh
import glob

datadir='../FTS data/15-02-5 150GHz Xpol/'
etalonlen = 3.25*25.4

def scan_analysis(dirpath, etalonlen, frequency, skiprows):

    myscan = fs.FTSscan(dirpath, frequency, etalonlen,skiprows=skiprows,\
        generateplots=True, useSincFitting=True, numinterppoints=15)
    myscan.initialize()
    myscan.driftcorrect()
    myscan.peakcorrect()
    myscan.symmetrize()
    myscan.getFFTs()
    myscan.averageFFT()
    myscan.getratio()
    myscan.checkguesses()
    myscan.fitparams()
    myscan.obtainerrorbars()
    myscan.savedata()

if __name__=='__main__':
    topdir = '../FTS data/'
    plotdir = 'plots/'
    scandirs = ['15-02-5 95GHz Ypol']
    skiprows = [24]
    freqs = [95]
    etalonlen = 3.25*25.4

    for i, scandir in enumerate(scandirs):
        print ("\nWorking on the scans {0}\n".format(scandir))
        dirpath = op.join(topdir, scandir, '')
        plotpath = op.join(plotdir,scandir, '')
        scan_analysis(dirpath, etalonlen, freqs[i], skiprows[i])

        print("\nAnalysis done. Moving all the plots to {0}\n".format(plotpath))
        if not op.isdir(plotpath):
            os.mkdir(plotpath)
        allplots = glob.glob('*.png')
        for plot in allplots:
            sh.mv("-f", plot, plotpath)

        