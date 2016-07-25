import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import glob
from astropy.constants import c
from scipy.fftpack import fft, fftfreq
c = c.value*1e-6 #Speed of light in mm/ns
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import matplotlib.mlab as mlab
import scipy.optimize as opt
import os.path as op
from scipy.stats import chi2
import emcee
import corner
import h5py
import scantypeenumerator as scantyper

# Fit functions for the interferogram peak
def quadraticfit(x, A, B, C):
    return A*x**2 + B*x + C

def sincfit(x, A, B, C, D):
    return A*(np.sinc(B*x+C))*np.cos(D*x)

def envelope(x, A, B, C):
    return A*(np.sinc(B*x+C)) 

def errorcheck(**kwargs):
    if (type(kwargs['datadir']) is not str):
        raise TypeError("The argument datadir must be a string\
         pointing to the location of the FTS Scans.")
    if (type(kwargs["frequency"]) is not int and type(kwargs["frequency"]) is not float):
        raise TypeError("The argument frequency must be a numerical value.")
    if (type(kwargs['etalonlength']) is not int and type(kwargs['etalonlength']) is not float):
        raise TypeError("The argument frequency must be a numerical value.")
    if (type(kwargs['skiprows']) is not int):
        raise TypeError("The argument skiprows must be an integer value.")
    if (type(kwargs['useonearm']) is not bool):
        raise TypeError("The argument useonearm must be True or False")
    if (type(kwargs['generateplots']) is not bool):
        raise TypeError("The argument generateplots must be True or False")
    if (type(kwargs['useSincFitting']) is not bool):
        raise TypeError("The argument useSincFitting must be True or False")

def get_save_names(dirname):
    save_name_list = op.basename(op.dirname(dirname)).split(' ')
    return '_'.join(save_name_list), '/' + '/'.join(save_name_list)

def loaddataset(self,filelist):
    ssignal = []
    sencoder = []
    for i, f in enumerate(filelist):
        index, time, encoder, signal = np.loadtxt(f, comments='#',\
         skiprows=self.skiprows, unpack=True)
        ssignal += [signal]
        sencoder += [encoder]
    print("All {0} files have been loaded".format(len(filelist)))
    ssignal = np.array(ssignal)
    sencoder = np.array(sencoder) 
    return ssignal, sencoder

def makeplots(self, Xlist, Ylist, tag, **kwargs):
    N = len(Xlist)
    for i in xrange(N):
        fig, ax = plt.subplots(figsize=(15,10))
        # print len(Xlist[i]), len(Ylist[i])    
        ax.plot(Xlist[i], Ylist[i])
        ax.set_xlabel(kwargs['x-label'])
        ax.set_ylabel(kwargs['y-label'])
        ax.axis('tight')
        plt.savefig(self.plt_savename +'-' + tag + '-' + kwargs['plt-type'] + str(i) + '.png')
        plt.close()

def convolutioncorrection(datascans, Nmask):
    signaldriftcorrected = []
    encoderdriftcorrected = []
    N = len(datascans['signal'])
    mask = np.ones(Nmask)
    startindex = Nmask-1
    endindex = -(Nmask-1)
    # indices = np.arange((datascans['signal'].shape[1]))
    # thresh = np.logical_and(indices >= Nmask-1, indices <= N - (Nmask-1))
    for i in xrange(N):
        signal = datascans['signal'][i]
        conv = np.convolve(datascans['signal'][i], mask, 'full')/Nmask
        # Region of full overlap where boundary effects are not visible
        signal = signal[startindex/2:endindex/2] - conv[startindex:endindex]
        signal -= np.average(signal)
        signaldriftcorrected += [signal]
        encoderdriftcorrected += [datascans['encoder'][i][startindex/2:endindex/2]]
    return np.array(signaldriftcorrected), np.array(encoderdriftcorrected)

def onearmcorrection(datascans, onearmscans):
    signaldriftcorrected = []
    encoderdriftcorrected = []
    N = len(datascans['signal'])
    for i in xrange(N):
        signal = (datascans['signal'][i] - onearmscans['signal-averaged'])
        signal -= np.average(signal)
        signaldriftcorrected += [signal]
    return np.array(signaldriftcorrected)

class FTSscan(object):
    def __init__(self, datadir, frequency, etalonlength, skiprows=0,\
     useonearm=True, generateplots=False,\
     useSincFitting=True):
        """

        """
        params = locals()
        errorcheck(**params)
        # Populate the class with its internal variables
        if (type(frequency) is float):
            frequency = int(frequency)
        self.fit95 = True if (frequency == 95) else False
        self.datadir = datadir
        self.frequency = frequency
        self.plt_savename, self.hdf5_name = get_save_names(datadir)
        self.useonearm = useonearm
        self.skiprows = skiprows
        self.etalonlength = etalonlength
        self.generateplots = generateplots
        # A series of flags to ensure pipeline is called in the right order
        self.initialized = False
        self.driftcorrected = False
        self.peakcorrected = False

        # I've written a routine that can traverse the directory and 
        # identify the relevant files that we want
        if (self.useonearm):
            self.sampleonearmls = scantyper.getsampleonearmscans(datadir)
            self.nosampleonearmls = scantyper.getnosampleonearmscans(datadir)
        else:
            self.sampleonearmls = []
            self.nosampleonearmls = []
        self.samplesls = scantyper.getsamplescans(datadir)
        self.nosamplesls = scantyper.getnosamplescans(datadir)

        self.filelist = glob.glob(self.datadir + '*.txt')
        if (self.filelist is []):
            raise IOError('file list of scans could not be generated')

    def initialize(self):
        signal, encoder = loaddataset(self, self.filelist)
        self.onearmdata = {}
        self.onearmdata['no-sample'] = {}
        self.onearmdata['sample'] = {}
        self.sampledata = {}
        self.nosampledata = {}

        if self.useonearm:
            self.onearmdata['no-sample']['signal'] = signal[self.nosampleonearmls]
            self.onearmdata['no-sample']['encoder'] = encoder[self.nosampleonearmls]
            self.onearmdata['sample']['signal'] = signal[self.sampleonearmls]
            self.onearmdata['sample']['encoder'] = encoder[self.sampleonearmls]
            
            self.onearmdata['no-sample']['signal-averaged']\
             = np.average(self.onearmdata['no-sample']['signal'], axis=0)
            self.onearmdata['no-sample']['encoder-averaged']\
             = np.average(self.onearmdata['no-sample']['encoder'], axis=0)
            self.onearmdata['sample']['signal-averaged']\
             = np.average(self.onearmdata['sample']['signal'], axis=0)
            self.onearmdata['sample']['encoder-averaged']\
             = np.average(self.onearmdata['sample']['encoder'], axis=0)
        self.nosampledata['signal'] = signal[self.nosamplesls]
        self.nosampledata['encoder'] = encoder[self.nosamplesls]
        self.sampledata['signal'] = signal[self.samplesls]
        self.sampledata['encoder'] = encoder[self.samplesls]
        # If generateplots is true, then generate plots here
        if self.generateplots:
            pltparams = {'x-label':r'Encoder [mm]',\
             'y-label':r'', 'plt-type':'raw-interferogram' }
            print ("Starting to make plots of raw interferogram ")
            makeplots(self, self.nosampledata['encoder'],\
             self.nosampledata['signal'], tag='no-sample', **pltparams)
            makeplots(self, self.sampledata['encoder'],\
             self.sampledata['signal'], tag='sample', **pltparams)
            if self.useonearm:
                makeplots(self, self.onearmdata['no-sample']['encoder'],\
                    self.onearmdata['no-sample']['signal'], tag='no-sample-one-arm', **pltparams)
                makeplots(self, self.onearmdata['sample']['encoder'],\
                    self.onearmdata['sample']['signal'], tag='sample-one-arm', **pltparams)
            print ("All the plots of the raw data have been completed ")
        print ("Scan successfully initialized")
        self.initialized = True

    def driftcorrect(self):
        if not self.initialized:
            raise RuntimeError("Scans must be initialized before drift correction")
        if self.useonearm:
            self.sampledata['signal-driftcorrected']\
             = onearmcorrection(self.sampledata, self.onearmdata['sample'])
            self.nosampledata['signal-driftcorrected']\
             = onearmcorrection(self.nosampledata, self.onearmdata['no-sample'])
            self.sampledata['encoder-driftcorrected'] = self.sampledata['encoder']
            self.nosampledata['encoder-driftcorrected'] = self.nosampledata['encoder']
        else:
            Nmask = 401
            self.sampledata['signal-driftcorrected'],\
             self.sampledata['encoder-driftcorrected'] = convolutioncorrection(self.sampledata, Nmask)
            self.nosampledata['signal-driftcorrected'],\
             self.nosampledata['encoder-driftcorrected'] = convolutioncorrection(self.nosampledata, Nmask)
        if self.generateplots:
            pltparams = {'x-label':r'Encoder [mm]',\
             'y-label':r'', 'plt-type':'driftcorrected-interferogram' }
            print ("Starting to make plots of the drift corrected interferograms")
            makeplots(self, self.nosampledata['encoder-driftcorrected'],\
             self.nosampledata['signal-driftcorrected'], tag='no-sample', **pltparams)
            makeplots(self, self.sampledata['encoder-driftcorrected'],\
             self.sampledata['signal-driftcorrected'], tag='sample', **pltparams)
            print ("All the plots of the raw data have been completed ")
        self.driftcorrected = True
        print ("Drift in source successfully corrected for")

    def peakcorrect(self):
        if not self.driftcorrected:
            raise RuntimeError("Scans must be drift corrected before peak correction")
        if self.useSincFitting:
            self.sampledata['signal-peakcorrected'],\
            self.sampledata['encoder-peakcorrected'] = sinccorrection(self.sampledata)
            self.nosampledata['signal-peakcorrected'],\
            self.nosampledata['encoder-peakcorrected'] = sinccorrection(self.nosampledata)
        else:
            self.sampledata['signal-peakcorrected'],\
            self.sampledata['encoder-peakcorrected'] = quadcorrection(self.sampledata)
            self.nosampledata['signal-peakcorrected'],\
            self.nosampledata['encoder-peakcorrected'] = quadcorrection(self.nosampledata)




