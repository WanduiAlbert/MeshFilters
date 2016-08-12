#!/usr/bin/env python

import numpy as np
import FTSscan as fs
import sys
import os.path as op
import os
import sh
import glob

def scan_analysis(dirpath, etalonlen, frequency):

    myscan = fs.FTSscan(dirpath, frequency, etalonlen,useonearm=False,\
        generateplots=True, useSincFitting=True, numinterppoints=15)
    myscan.initialize()
    myscan.driftcorrect()
    myscan.peakcorrect()
    myscan.symmetrize()
    # myscan.generateplots = True 
    myscan.getFFTs()
    # myscan.generateplots = False 
    myscan.averageFFT()
    myscan.getratio()
    # myscan.generateplots = False
    myscan.checkguesses()
    # myscan.generateplots = True 
    myscan.fitparams()
    myscan.obtainerrorbars()
    myscan.savedata()

if __name__=='__main__':
    topdir = '../FTS data/'
    plotdir = 'plots/'
    scandirs = ['18-01 150GHz Xpol', '18-01 150GHz Ypol',\
    '18-01 95GHz Xpol', '18-01 95GHz Ypol']
    freqs = [150, 150, 95, 95]
    etalonlen = 2.25*25.4#2.625*25.4

    for i, scandir in enumerate(scandirs):
        print ("\nWorking on the scans {0}\n".format(scandir))
        dirpath = op.join(topdir, scandir, '')
        plotpath = op.join(plotdir,scandir, '')
        scan_analysis(dirpath, etalonlen, freqs[i])

        print("\nAnalysis done. Moving all the plots to {0}\n".format(plotpath))
        if not op.isdir(plotpath):
            os.mkdir(plotpath)
        sh.mv("-f", sh.glob('*.png'), plotpath)

        