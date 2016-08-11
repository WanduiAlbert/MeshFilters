import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import glob
from astropy.constants import c
from scipy.fftpack import fft, fftfreq
c = c.value*1e-6 #Speed of light in mm/ns
from scipy.interpolate import splev, splrep
from scipy.optimize import curve_fit
import matplotlib.mlab as mlab
import scipy.optimize as opt
import os.path as op
from scipy.stats import chi2
import emcee
import corner
import h5py
import scantypeenumerator as scantyper
from scipy.integrate import simps

elength = 0
# Fit functions for the interferogram peak
def quadraticfit(x, A, B, C):
    return A*x**2 + B*x + C

def sincfit(x, A, B, C, D):
    return A*(np.sinc(B*x+C))*np.cos(D*x)

def envelope(x, A, B, C):
    return A*(np.sinc(B*x+C)) 

# Functions for fitting for the parameters
def transmissionModel(nu, R, T, L, C):
    A = (T/(1-R))**2
    F = 4*R/(1-R)**2
    B = 2*2*np.pi/c*L*elength
    return A*1./(1 + F * np.sin((B*nu + C*np.pi)/2)**2)

def lnlike(theta, x, y, yerr):
    model = transmissionModel(x, *theta)
    inv_sigma2 = 1.0/(yerr**2)
    return -0.5*(np.sum((y-model)**2*inv_sigma2))

def chisq(theta, x, y, yerr):
    N = len(x)
    return -2.0*lnlike(theta, x, y, yerr)

def pvalue(chisqval, dof):
    return chi2.sf(chisqval, dof)

nll = lambda *args: -lnlike(*args) #Negative Log Likelihood

def lnprior(theta):
    R, T, L, C = theta
    if (0.0 < R < 1.0 and
        0.0 < T < 1.0 and
        # 0.0 < R + T <= 1.0 and 
        0.5 < L < 1.5 and
        -2.0 < C < 2.0) :
        return 0.0
    else:
        return -np.inf

def lnprob(theta, x, y, yerr):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, yerr)

def errorcheck(**kwargs):
    if (type(kwargs['datadir']) is not str):
        raise TypeError("The argument datadir must be a string\
         pointing to the location of the FTS Scans.")
    if (type(kwargs["frequency"]) is not int and type(kwargs["frequency"]) is not float):
        raise TypeError("The argument frequency must be a numerical value.")
    if (type(kwargs['etalonlength']) is not int and type(kwargs['etalonlength']) is not float):
        raise TypeError("The argument frequency must be a numerical value.")
    if (type(kwargs['useonearm']) is not bool):
        raise TypeError("The argument useonearm must be True or False")
    if (type(kwargs['generateplots']) is not bool):
        raise TypeError("The argument generateplots must be True or False")
    if (type(kwargs['useSincFitting']) is not bool):
        raise TypeError("The argument useSincFitting must be True or False")

def get_save_names(dirname):
    save_name_list = op.basename(op.dirname(dirname)).split(' ')
    return '_'.join(save_name_list), '/' + '/'.join(save_name_list)

def skipstart(filepath):
    f = open(filepath, 'r')
    while True:
        line = f.readline().strip()
        if (line == '# DATA'):
            break
    return f

def loaddataset(self,filelist):
    ssignal = []
    sencoder = []
    for i, filepath in enumerate(filelist):
        f = skipstart(filepath)
        index, time, encoder, signal = np.loadtxt(f, comments='#',\
            unpack=True)
        ssignal += [signal]
        sencoder += [encoder]
    print("All {0} files have been loaded".format(len(filelist)))
    ssignal = np.vstack(ssignal)
    sencoder = np.vstack(sencoder) 
    return ssignal, sencoder

def makeplots(self, Xlist, Ylist, tag, **kwargs):
    N = len(Ylist)
    for i in xrange(N):
        fig, ax = plt.subplots(figsize=(15,10))
        # print len(Xlist[i]), len(Ylist[i]) 
        ax.plot(Xlist[i], Ylist[i])
        ax.set_xlabel(kwargs['x-label'])
        ax.set_ylabel(kwargs['y-label'])
        ax.axis('tight')
        ax.set_xticklabels(["{0:3.1f}".format(t) for t in ax.get_xticks()])
        ax.set_yticklabels(["{0:1.6f}".format(t) for t in ax.get_yticks()])
        plt.savefig(self.plt_savename +'-' + tag + '-' + kwargs['plt-type'] + str(i) + '.png')
        plt.close()

def fftplot(self, datascans, tag, **kwargs):
    fig, ax = plt.subplots(figsize=(15,10))
    N = len(datascans['signal'])
    for i in xrange(N):
        ax.plot(datascans['freq-interest'][i][::2], datascans['fft-interest'][i][::2], label=str(i))
        ax.set_xlabel(kwargs['x-label'])
        ax.set_ylabel(kwargs['y-label'])
        ax.axis('tight')
    ax.errorbar(self.frequency[self.thresh][::2], datascans['fft-averaged'][::2],\
        yerr = datascans['fft-error'][::2], label="averaged",fmt='-',\
        ecolor='k', color='k')
    ax.legend(loc='best')
    plt.savefig(self.plt_savename +'-' + tag + '-' + kwargs['plt-type'] + '.png')
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
        conv = np.convolve(signal, mask, 'full')/Nmask
        # Region of full overlap where boundary effects are not visible
        signal = signal[startindex/2:endindex/2] - conv[startindex:endindex]
        signal -= np.average(signal)
        signaldriftcorrected += [signal]
        encoderdriftcorrected += [datascans['encoder'][i][startindex/2:endindex/2]]
    return np.array(signaldriftcorrected), np.array(encoderdriftcorrected)

def onearmcorrection(datascans, onearmscans):
    signaldriftcorrected = []
    N = len(datascans['signal'])
    for i in xrange(N):
        signal = (datascans['signal'][i] - onearmscans['signal-averaged'])
        # signal /= np.sqrt(onearmscans['signal-averaged'])
        signal -= np.average(signal)
        signaldriftcorrected += [signal]
    return np.array(signaldriftcorrected)

def correction(datascans, fitfunction, initialparams, minindex, maxindex):
    popts = []
    pcovs = []
    N = len(datascans['signal'])
    for i in xrange(N):
        x = datascans['encoder-driftcorrected'][i][minindex[i]:maxindex[i]]
        y = datascans['signal-driftcorrected'][i][minindex[i]:maxindex[i]]
        popt, pcov = curve_fit(fitfunction,x,y,p0=initialparams)
        xnew = np.r_[250:550:1000j]
        ynew= envelope(xnew,*popt[:-1])
        yguess = envelope(xnew, *initialparams[:-1])
        fig, ax = plt.subplots(figsize=(15,10))
        ax.plot(datascans['encoder-driftcorrected'][i], datascans['signal-driftcorrected'][i], 'b.-')
        ax.plot(xnew, ynew, 'r-')
        ax.plot(xnew, yguess, 'k-')
        ax.grid(which='major', axis='x', linewidth=0.75, linestyle='-', color='0.95')
        ax.grid(which='minor', axis='x', linewidth=0.25, linestyle='-', color='0.95')
        ax.grid(which='major', axis='y', linewidth=0.75, linestyle='-', color='0.95')
        ax.grid(which='minor', axis='y', linewidth=0.25, linestyle='-', color='0.95')
        ax.axis('tight')
        plt.savefig(str(i) + '.png')
        plt.close()
        pcovs += [pcov]
        popts += [popt]
    return popts, pcovs

def sinccorrection(self, datascans):
    initialparams = []
    if self.fit95:
        initialparams = np.array([1.5e-4,40/c,-51.6,1190./c]) # When no sqrt
        # initialparams = np.array([1.4e-2,40/c,-51.6,1190./c])
    else:
        initialparams = np.array([1.5e-4,40/c,-51.6,1862.3/c]) # When no sqrt
        # initialparams = np.array([1.4e-2,40/c,-51.6,1862.3/c])
    minindex = map(lambda x: np.where(x >= 375)[0][0], datascans['encoder-driftcorrected'])
    maxindex = map(lambda x: np.where(x <= 400)[0][-1], datascans['encoder-driftcorrected'])
    # print(minindex, maxindex)
    # print (datascans['encoder'][0][minindex[0]])

    popts, pcovs = correction(datascans, sincfit, initialparams, minindex, maxindex)
    signalpeaks = np.array(map(lambda x: sincfit(-x[2]/x[1], *x), popts))
    cosineshift = np.array(map(lambda x: np.cos((-x[2]/x[1])*x[-1]), popts))
    phaseshift = np.round(map(lambda x: ((-x[2]/x[1])*x[-1])/(np.pi), popts))
    encoderpeaks = np.array(map(lambda x, N: np.pi/x[-1]*N, popts, phaseshift))[:,np.newaxis]

    print (encoderpeaks.shape)

    self.signalpeaks = signalpeaks
    # encoderpeaks = np.array(map(lambda x:-x[2]/x[1], popts))[:,np.newaxis]
    relerr = np.array(map(lambda sigx, x: (sigx[2,2]/(x[2]**2) +
        sigx[1,1]/(x[1])**2 - 2*sigx[2,1]/(x[1]*x[2]))**0.5, pcovs, popts))[:,np.newaxis]
    encoderpeakerr = encoderpeaks*relerr
    print("\n")
    for i in xrange(len(encoderpeaks)):
        # print ("The cos at zero pd is {0:1.5f} with a phase shift {1:1.5f}".format(cosineshift[i], phaseshift[i]))
        print("for the scan {0:d} the position ".format(i)+\
            "of zero p.d is {0:1.6f} +/- {1:1.6f}".format(encoderpeaks[i,0], encoderpeakerr[i,0]))

    # for i in xrange(len(encoderpeaks)):
    #     xnew = np.r_[375:400:1000j]
    #     ynew= sincfit(xnew,*popts[i])
    #     # yguess = envelope(xnew, *initialparams[:-1])
    #     fig, ax = plt.subplots(figsize=(15,10))
    #     ax.plot(datascans['encoder-driftcorrected'][i][minindex[i]:maxindex[i]] - encoderpeaks[i],\
    #         datascans['signal-driftcorrected'][i][minindex[i]:maxindex[i]], 'b.-')
    #     ax.plot(xnew - encoderpeaks[i], ynew, 'r-')
    #     # ax.plot(xnew, yguess, 'k-')
    #     ax.grid(which='major', axis='x', linewidth=0.75, linestyle='-', color='0.00')
    #     ax.grid(which='minor', axis='x', linewidth=0.25, linestyle='-', color='0.00')
    #     ax.grid(which='major', axis='y', linewidth=0.75, linestyle='-', color='0.00')
    #     ax.grid(which='minor', axis='y', linewidth=0.25, linestyle='-', color='0.00')
    #     ax.axis('tight')
    #     plt.savefig(str(i) + '.png')
    #     plt.close()
    return (datascans['encoder-driftcorrected'] - encoderpeaks)*(2/c)

def quadcorrection(datascans):
    peaks = np.argmax(datascans['signal-driftcorrected'], axis=1)
    minindex = peaks - 1
    maxindex = peaks + 2
    popts, pcovs = correction(datascans, quadraticfit, [1,1,1], minindex, maxindex)
    encoderpeaks = np.array(map(lambda x:-x[1]/(2*x[0]), popts))[:,np.newaxis]
    relerr = np.array(map(lambda sigx, x: ((sigx[1]/x[1])**2 + (sigx[0]/x[0])**2\
        - 2*sigx[1,0]/(x[1]*x[0]))**0.5, pcovs, popts))[:,np.newaxis]
    encoderpeakerr = encoderpeaks*relerr
    for i in xrange(len(encoderpeaks)):
        print("for the scan {0:d} the position\
            of zero p.d is {1:1.6f} +/- {2:1.6f}".format(i, encoderpeaks[i], encoderpeakerr[i]))
    return (datascans['encoder-driftcorrected'] - encoderpeaks)*(2/c)

def resamplesig(datascans, x_new):
    signew=[]
    N = len(datascans['signal'])

    for i in xrange(N):
        tck, fp, ier, msg = splrep(datascans['encoder-driftcorrected'][i],\
         datascans['signal-driftcorrected'][i], w=None,k=3, task=0,s=None,full_output=True)
        print (msg)
        y = splev(x_new, tck, der=0, ext=1)
        signew += [y]
        print( "Finished working on the %d scan" %i)
    return np.array(signew)

def getfft(datascans):
    N = len(datascans['signal'])
    ffts = []
    for i in xrange(N):
        y = fft((datascans['signal-resampled'][i]))
        ffts += [y] 
    return np.array(ffts)

def getthresh(self, datascans):
    N = len(datascans['signal'])
    x = np.vstack([self.frequency[self.thresh]]*N)
    y = (np.real(datascans['signal-fft'])[np.vstack([self.thresh]*N)])
    y = y.reshape((N,-1))
    return x, y

def validateguess(rawguess):
    try:
        guesses = [float(i) for i in rawguess]
        print ("The new guesses are {0}".format(guesses))
        return guesses
    except ValueError:
        print ('The guesses supplied must be numerical values')
        return None

def getpower(y, x):
    return simps(y**2, x=x, even='avg', axis=1)

def obtainguess():
    guesses = [

    ]
    while (True):
        try:
            rawguess = raw_input("\nEnter an updated guess in the format R, T, L, C: ").split(',')
            guesses = validateguess(rawguess)
            if guesses is not None: break
        except EOFError:
            print ('Proceeding with current guess')
            break
    return guesses

def recursivesave(hdf5file, mydict):
    for key, item in mydict.items():
        if isinstance(item, dict):
            try:
                curr_grp = hdf5file[key]
            except (KeyError, RuntimeError):
                curr_grp = hdf5file.create_group(key)
            recursivesave(curr_grp, item)
        else:
            savedataset(hdf5file, item, key)

def savedataset(hdf5file, data, key):
    data = np.array(data)
    try:
        hdf5file[key][...] = data
    except (KeyError, RuntimeError):
        hdf5file.create_dataset(key, data= data, maxshape=(100, 10000))
    # except TypeError:
    #     del hdf5file[key]
    #     hdf5file.create_dataset(key, data= data, maxshape=(100, 10000))

class FTSscan(object):
    def __init__(self, datadir, frequency, etalonlength,\
        useonearm=True, generateplots=False,\
        useSincFitting=True, numinterppoints=15):
        """

        """
        params = locals()
        errorcheck(**params)
        # Populate the class with its internal variables
        if (type(frequency) is float):
            frequency = int(frequency)
        self.fit95 = True if (frequency == 95) else False
        self.datadir = op.join(datadir, '')
        self.frequency = frequency
        self.plt_savename, self.hdf5_name = get_save_names(datadir)
        self.useonearm = useonearm
        global elength
        self.etalonlength = etalonlength
        elength = self.etalonlength # Wierd global!!!
        self.generateplots = generateplots
        self.useSincFitting = useSincFitting
        self.numinterppoints = numinterppoints
        self.dx = 0.0

        # A series of flags to ensure pipeline is called in the right order
        self.initialized = False
        self.driftcorrected = False
        self.peakcorrected = False
        self.symmetrized = False
        self.transformed = False
        self.averaged = False
        self.ratioed = False
        self.guesschecked = False
        self.fitted = False

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

        print ("SAMPLE - ONE ARM SCANS: {0}".format(self.sampleonearmls))
        print ("NO SAMPLE - ONE ARM SCANS: {0}".format(self.nosampleonearmls))
        print ("SAMPLE: {0}".format(self.samplesls))
        print ("NO SAMPLE: {0}".format(self.nosamplesls))

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
            self.onearmdata['no-sample']['signal'] = signal[self.sampleonearmls]
            self.onearmdata['no-sample']['encoder'] = encoder[self.sampleonearmls]
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
            print ("All the plots of the drift corrected interferograms have been completed ")
        self.driftcorrected = True

        print ("Drift in source successfully corrected for")

    def peakcorrect(self):
        if not self.driftcorrected:
            raise RuntimeError("Scans must be drift corrected before peak correction")
        if self.useSincFitting:
            self.sampledata['encoder-driftcorrected'] = sinccorrection(self, self.sampledata)
            self.nosampledata['encoder-driftcorrected'] = sinccorrection(self, self.nosampledata)
        else:
            self.sampledata['encoder-driftcorrected'] = quadcorrection(self.sampledata)
            self.nosampledata['encoder-driftcorrected'] = quadcorrection(self.nosampledata)

        print ("\nChecking for the power in the sample vs no sample:\n")
        y = self.sampledata['signal-driftcorrected']
        x = self.sampledata['encoder-driftcorrected']
        print (x.shape)
        print (y.shape)
        self.sampledata['power'] = getpower(y,x)

        y = self.nosampledata['signal-driftcorrected']
        x = self.nosampledata['encoder-driftcorrected']
        self.nosampledata['power'] = getpower(y,x)

        print ("For sample: Average Power is {0}\n".format(np.average(self.sampledata['power'])))
        print ("For no sample: Average Power is {0}\n".format(np.average(self.nosampledata['power'])))

        if self.generateplots:
            pltparams = {'x-label':r'Path Difference [ns]',\
             'y-label':r'', 'plt-type':'peakcorrected-interferogram' }
            print ("Starting to make plots of the peak corrected interferograms")
            makeplots(self, self.nosampledata['encoder-driftcorrected'],\
             self.nosampledata['signal-driftcorrected'], tag='no-sample', **pltparams)
            makeplots(self, self.sampledata['encoder-driftcorrected'],\
             self.sampledata['signal-driftcorrected'], tag='sample', **pltparams)
            print ("All the plots of the peak corrected interferograms have been completed ")
        self.peakcorrected= True

    def symmetrize(self):
        if not self.peakcorrected:
            raise RuntimeError("Scans must be peak corrected before resampling")
        n = self.numinterppoints
        k = 10.0
        self.dx  = k*2**(1 - n)
        encoder_new = self.dx * np.arange(0, 2**n + 1) - k
        self.encoder_resampled= encoder_new[:-1]
        self.sampledata["signal-resampled"] = resamplesig(self.sampledata, self.encoder_resampled)
        self.nosampledata["signal-resampled"] = resamplesig(self.nosampledata, self.encoder_resampled)
        geteven = lambda x: 0.5*(x[1:] + x[::-1][:-1])
        getodd = lambda x: 0.5*(x[1:] - x[::-1][:-1])

        # # Print out the value at the peak both before and after resampling
        # for i in xrange(len(self.sampledata['signal'])):
        #     print ("\nBefore resampling, the peak value was {0:1.9f}".format(self.signalpeaks[i]))
        #     print ("After resampling,"+\
        #         " the peak value was {0:1.9f}\n".format(self.sampledata['signal-resampled'][i][2**(n-1)]))

        if self.generateplots:
            fig, ax = plt.subplots(figsize=(15,10))
            ax.plot(self.encoder_resampled[1:],\
                geteven(self.nosampledata['signal-resampled'][0]),\
                'r-', label='Even')
            ax.plot(self.encoder_resampled[1:],\
                getodd(self.nosampledata['signal-resampled'][0]),\
                'b-', label='Odd')
            ax.legend(loc='best')
            plt.savefig('nosampledata_evenodd.png')

            fig, ax = plt.subplots(figsize=(15,10))
            ax.plot(self.encoder_resampled[1:],\
                geteven(self.sampledata['signal-resampled'][0]),\
                'r-', label='Even')
            ax.plot(self.encoder_resampled[1:],\
                getodd(self.sampledata['signal-resampled'][0]),\
                'b-', label='Odd')
            ax.legend(loc='best')
            plt.savefig('sampledata_evenodd.png')

            # pltparams = {'x-label':r'Path Difference [ns]',\
            #  'y-label':r'Signal', 'plt-type':'resampled-interferogram' }
            # print ("Starting to make plots of the resampled interferograms")
            # makeplots(self, [self.encoder_resampled]*len(self.nosampledata['signal']),\
            #  self.nosampledata['signal-resampled'], tag='no-sample', **pltparams)
            # makeplots(self, [self.encoder_resampled]*len(self.sampledata['signal']),\
            #  self.sampledata['signal-resampled'], tag='sample', **pltparams)
            print ("All the plots of the symmetrized interferograms have been completed ")
        self.symmetrized = True

    def getFFTs(self):
        if not self.symmetrized:
            raise RuntimeError("Scans must be symmetrized before taking FFTs")
        self.frequency = fftfreq(len(self.encoder_resampled), self.dx)
        self.sampledata['signal-fft'] = getfft(self.sampledata)
        self.nosampledata['signal-fft'] = getfft(self.nosampledata)
        if self.fit95:
            self.thresh = np.logical_and(self.frequency >= 90,\
             self.frequency <= 100) #For 95 GHz
        else:
            self.thresh = np.logical_and(self.frequency >= 140,\
             self.frequency <= 159)

        #sample case
        self.nosampledata['freq-interest'], self.nosampledata['fft-interest'] = getthresh(self, self.nosampledata)
        self.sampledata['freq-interest'], self.sampledata['fft-interest'] = getthresh(self, self.sampledata)

        x = self.frequency
        self.sampledata['fft-power'] = getpower(np.real(self.sampledata['signal-fft']), x)
        self.nosampledata['fft-power'] = getpower(np.real(self.nosampledata['signal-fft']), x)

        x = self.frequency[self.thresh]
        self.sampledata['fft-power-interest'] = getpower(np.real(self.sampledata['fft-interest']), x)
        self.nosampledata['fft-power-interest'] = getpower(np.real(self.nosampledata['fft-interest']), x)
        
        print ("For sample: Average Power is {0} over all frequencies.\n".\
            format(np.average(self.sampledata['fft-power'])))
        print ("For no sample: Average Power is {0} over all frequencies.\n".\
            format(np.average(self.nosampledata['fft-power'])))
        print ("For sample: Average Power is {0} over our freq band.\n".\
            format(np.average(self.sampledata['fft-power-interest'])))
        print ("For no sample: Average Power is {0} over our freq band.\n".\
            format(np.average(self.nosampledata['fft-power-interest'])))

        if self.generateplots:
            
            pltparams = {'x-label':r'Frequency [GHz]',\
             'y-label':r'Signal', 'plt-type':'fft' }
            
            print ("Starting to make plots of the fourier transforms")
            # fig, ax = plt.subplots(figsize=(15,10))
            # ax.plot(self.sampledata['freq-interest'][0][::4],\
            #     np.real(self.sampledata['signal-fft'][0][self.thresh][::4]),\
            #     'r-', label='Real')
            # ax.plot(self.sampledata['freq-interest'][0][::4],\
            #     np.abs(np.imag(self.sampledata['signal-fft'][0][self.thresh][::4])),\
            #     'b-', label='Imag')
            # ax.legend(loc='best')
            # plt.savefig('sampledata_realvsimag.png')

            # fig, ax = plt.subplots(figsize=(15,10))
            # ax.plot(self.nosampledata['freq-interest'][0][::4],\
            #     np.real(self.nosampledata['signal-fft'][0][self.thresh][::4]),\
            #     'r-', label='Real')
            # ax.plot(self.nosampledata['freq-interest'][0][::4],\
            #     np.abs(np.imag(self.nosampledata['signal-fft'][0][self.thresh][::4])),\
            #     'b-', label='Imag')
            # ax.legend()
            # plt.savefig('nosampledata_realvsimag.png')
            # plt.close()
            makeplots(self, self.nosampledata['freq-interest'],\
             self.nosampledata['fft-interest'],tag='no-sample', **pltparams)

            makeplots(self, self.sampledata['freq-interest'],\
             self.sampledata['fft-interest'],tag='sample', **pltparams)
            print ("All the plots of the fourier transforms have been completed ")
        self.transformed = True

    def averageFFT(self):
        if not self.transformed:
            raise RuntimeError("FFTs of scans must be taken before they can be averaged")
        self.sampledata['fft-averaged'] = np.mean(np.real(self.sampledata['fft-interest']), axis=0)
        self.sampledata['fft-error'] = np.std(np.real(self.sampledata['fft-interest']), axis=0)
        self.nosampledata['fft-averaged'] = np.mean(np.real(self.nosampledata['fft-interest']), axis=0)
        self.nosampledata['fft-error'] = np.std(np.real(self.nosampledata['fft-interest']), axis=0)

        if self.generateplots:
            pltparams = {'x-label':r'Frequency [GHz]',\
             'y-label':r'Signal', 'plt-type':'averaged-fft' }

            print ("\nStarting to make plots of the fourier transforms")
            fftplot(self, self.nosampledata,tag='no-sample', **pltparams)

            fftplot(self, self.sampledata,tag='sample', **pltparams)
            
            print ("All the plots of the FFT have been completed")

            # Plot the averaged sample vs no-sample
            fig, ax = plt.subplots(figsize=(15,10))
            ax.plot(self.frequency[self.thresh][::2], self.nosampledata['fft-averaged'][::2], 'k', label='no sample')
            ax.plot(self.frequency[self.thresh][::2], self.sampledata['fft-averaged'][::2], 'r', label='sample');
            ax.set_xlabel('Frequency [GHz]')
            ax.set_ylabel('Signal');
            ax.legend(loc='best');
            ax.axis('tight');
            plt.savefig(self.plt_savename + 'no-sample_vs_sample_fft.png')
            plt.close()
        self.averaged = True

    def getratio(self):
        if not self.averaged:
            raise RuntimeError("FFTs of scans must be averaged before the ratio can be computed")
        self.nosampledata['relative-error'] = self.nosampledata['fft-error']/self.nosampledata['fft-averaged']
        self.sampledata['relative-error'] = self.sampledata['fft-error']/self.sampledata['fft-averaged']
        # self.ratio = self.sampledata['fft-interest']/self.nosampledata['fft-interest']
        self.ratio_avg= self.sampledata['fft-averaged']/self.nosampledata['fft-averaged']
        self.ratio_err = self.ratio_avg*np.sqrt(self.nosampledata['relative-error']**2 + \
            self.sampledata['relative-error']**2)

        # f = open('no-sqrt-1801.txt', 'w')
        # f.write('# Frequency Ratio Error\n')
        # np.savetxt(f, np.c_[self.frequency[self.thresh], self.ratio_avg, self.ratio_err], delimiter=' ')

        if self.generateplots:
            # Plot the relative error
            print("\nStarting to make plots for the transmission ratio")
            fig, ax = plt.subplots(figsize=(15,10))
            ax.plot(self.frequency[self.thresh][::2], self.sampledata['relative-error'][::2], 'b', label='Sample')
            ax.plot(self.frequency[self.thresh][::2], self.nosampledata['relative-error'][::2], 'g', label='No Sample')
            ax.legend(loc='best');
            ax.axis('tight');
            ax.set_xlabel('Frequency [GHz]')
            plt.savefig(self.plt_savename + 'relative_error.png')


            # Plot the ratios from all the plots
            pltparams = {'x-label':r'Frequency [GHz]',\
             'y-label':r'Signal', 'plt-type':'ratios' }
            fig, ax = plt.subplots(figsize=(15,10))
            # N = len(datascans['signal'])
            # for i in xrange(N):
            #     ax.plot(self.sampledata['freq-interest'], self.ratio, label=str(i))
            #     ax.set_xlabel(kwargs['x-label'])
            #     ax.set_ylabel(kwargs['y-label'])
            #     ax.axis('tight')
            ax.errorbar(self.frequency[self.thresh][::2], self.ratio_avg[::2],\
                yerr = self.ratio_err[::2],fmt='-',\
                ecolor='k', color='k')
            xmin, xmax, ymin, ymax = ax.axis('tight')
            ax.hlines(1.0, xmin, xmax, colors='r', linestyles='solid');
            plt.savefig(self.plt_savename +'-' + 'ratio-averaged' + '.png')
            plt.close()
            print ("All ratio plots completed")
        self.ratioed = True

    def checkguesses(self):
        if not self.ratioed:
            raise RuntimeError("The ratio must be computed before the guesses can be checked")
        fig, ax = plt.subplots(figsize=(15,10))
        pltparams = {'x-label':r'Frequency [GHz]',\
        'y-label':r'Signal', 'plt-type':'ratios' }
        R = 0.02
        L = np.cos(3*np.pi/180)
        C = 1.1
        T = 0.98

        guesses = [R, T, L, C]
        print ("\nStarting with the initial guesses {0}".format(guesses))
        while (True):
            print ("\nGenerating plot using guesses {0}".format(guesses))
            plt.cla()
            ax.errorbar(self.frequency[self.thresh][::2], self.ratio_avg[::2],\
                yerr = self.ratio_err[::2],fmt='-',\
                ecolor='b', color='b', label='data')
            ax.plot(self.frequency[self.thresh],\
             transmissionModel(self.frequency[self.thresh], *guesses),\
             'r-', label='Guessed Fit');
            ax.legend(loc='best')
            ax.axis('tight');
            ax.set_xticklabels(["{0:3.1f}".format(t) for t in ax.get_xticks()])
            ax.set_yticklabels(["{0:1.4f}".format(t) for t in ax.get_yticks()])
            ax.grid(which='major')
            ax.set_xlabel('Frequency [GHz]')
            plt.savefig(self.plt_savename +'-' + 'guessed-fit' + '.png')

            try: 
                response = raw_input("Is the guess of the parameters acceptable? ")
                if response.lower() in ['y', 'yes']:
                    break
                elif response.lower() in ['n', 'no']:
                    newguess = obtainguess()
                    if not newguess:
                        break
                    else:
                        guesses = newguess
                else:
                    print ('Enter yes or no')
            except EOFError:
                print ('Proceeding with current guess')  
                break
        
        # Make final plot with updated guesses
        plt.cla()
        ax.errorbar(self.frequency[self.thresh][::2], self.ratio_avg[::2],\
            yerr = self.ratio_err[::2],fmt='-',\
            ecolor='b', color='b', label='data')
        ax.plot(self.frequency[self.thresh],\
            transmissionModel(self.frequency[self.thresh], *guesses),\
            'r-', label='Guessed Fit');
        ax.legend(loc='best')
        ax.axis('tight');
        ax.grid(which='major')
        ax.set_xticklabels(["{0:3.1f}".format(t) for t in ax.get_xticks()])
        ax.set_yticklabels(["{0:1.4f}".format(t) for t in ax.get_yticks()])
        ax.set_xlabel('Frequency [GHz]')
        plt.savefig(self.plt_savename +'-' + 'guessed-fit' + '.png')
        plt.close()
        self.guesses = guesses
        self.guesschecked = True
     
    def fitparams(self):
        if not self.guesschecked:
            raise RuntimeError("Check the guesses of the parameters before attempting to fit")
        N = len(self.frequency[self.thresh])
        self.dof = N - len(self.guesses)

        result = opt.minimize(nll, self.guesses, args=(self.frequency[self.thresh],\
            self.ratio_avg,self.ratio_err), method='Powell')
        self.params = result['x']
        self.chisqval = chisq(self.params, self.frequency[self.thresh],\
                    self.ratio_avg,self.ratio_err)
        self.chisqreduced = self.chisqval/self.dof
        self.pval = pvalue(self.chisqval, self.dof) 

        resultstr = "\nParameters extracted with a reduced chi-square "
        resultstr += "value of {0:4.5f} and a p-value {1:4.5f}".format(self.chisqreduced, self.pval)
        print (result['message'] + resultstr)

        paramstr = " R = {params[0]}, T = {params[1]}, L= {params[2]}, C = {params[3]}".format(params=self.params)
        print "optimal parameters are:" + paramstr

        if self.generateplots:
            print ("Generating a plot of the best fit to the data")
            fig, ax = plt.subplots(figsize=(15,10))
            ax.errorbar(self.frequency[self.thresh][::2], self.ratio_avg[::2],\
                yerr = self.ratio_err[::2],fmt='-',\
                ecolor='b', color='b', label='data')
            ax.plot(self.frequency[self.thresh],\
             transmissionModel(self.frequency[self.thresh], *self.params),\
             'r-', label='Best Fit');
            ax.legend(loc='best')
            ax.axis('tight');
            ax.grid(which='major')
            ax.set_xlabel('Frequency [GHz]')
            ax.set_xticklabels(["{0:3.1f}".format(t) for t in ax.get_xticks()])
            ax.set_yticklabels(["{0:1.4f}".format(t) for t in ax.get_yticks()])
            plt.savefig(self.plt_savename +'-' + 'best-fit' + '.png')
            plt.close()
            print ("Fitting of the parameters completed. Moving to calculating the errorbars!")
        self.fitted = True

    def obtainerrorbars(self):
        if not self.fitted:
            raise RuntimeError("The best fit must be computed before the error on the parameters estimated")
        labels = [r'R', r'T', r'L', r'C']
        ndim, nwalkers = 4, 500
        pos = [self.params + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(self.frequency[self.thresh],\
            self.ratio_avg,self.ratio_err))
        nsteps = 2000
        sampler.run_mcmc(pos, nsteps);

        if self.generateplots:
            try:
                print ("\nGenerating plots of the movement of walkers")
                # Show the movement of the walkers in paramspace
                for j in xrange(ndim):
                    fig,ax = plt.subplots(figsize=(15,10))
                    for i in xrange(0, nwalkers, 4):
                        if i is 0:
                            ax.plot(sampler.chain[i,:,j], 'k', label='Simulated Fit')
                        ax.plot(sampler.chain[i,:,j], 'k')
                    __, __, ymin, ymax = ax.axis()
                    ax.vlines(200, ymin,ymax,color=u'r')
                    ax.axis('tight')
                    ax.grid(which='both');
                    ax.set_ylabel(labels[j])
                    ax.set_xlabel(r'Number of steps')
                    ax.set_xticklabels(["{0:3.1f}".format(t) for t in ax.get_xticks()])
                    ax.set_yticklabels(["{0:1.4f}".format(t) for t in ax.get_yticks()])
                    plt.savefig(self.plt_savename + '{0}'.format(labels[j]) + '-walkers-plot.png')
                    plt.close()
            except RuntimeError as err:
                print ("RuntimeError: {0}.\nPlot wasn't generated. Continuing".format(err))

        samples = sampler.chain[:, 200:, :].reshape((-1, ndim))

        if self.generateplots:
            print ("Generating a corner plot displaying the results")
            try:
                # Make a corner plot of the results
                figure, ax = plt.subplots(nrows = ndim, ncols=ndim, figsize=(15,12))
                fig = corner.corner(samples, title_fmt='.4f',labels=["$R$", "$T$", "$L$", "$C$"],\
                    truths=self.params, quantiles=[.16, .50, .84],\
                    fig=figure,show_titles=True, use_math_text=True)
                fig.savefig(self.plt_savename + 'corner_plot.png')
                plt.close()
            except RuntimeError as err:
                print ("RuntimeError: {0}.\nPlot wasn't generated. Continuing".format(err))

        self.R_mcmc, self.T_mcmc, self.L_mcmc, self.C_mcmc =\
        map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),\
            zip(*np.percentile(samples, [16, 50, 84],axis=0)))

        self.A_mcmc = (1-(self.R_mcmc[0] + self.T_mcmc[0]),\
            ((self.R_mcmc[1]**2 + self.T_mcmc[1]**2))**0.5,\
            (self.R_mcmc[2]**2 + self.T_mcmc[2]**2)**0.5)

        print ("\nWe determined R to be {0:1.5f} plus {1:1.5f} or minus {2:1.5f}".format(*self.R_mcmc))
        print ("We determined T to be {0:1.5f} plus {1:1.5f} or minus {2:1.5f}".format(*self.T_mcmc))
        print ("We determined A to be {0:1.5f} plus {1:1.5f} or minus {2:1.5f}".format(*self.A_mcmc))
        print ("We determined L to be {0:1.5f} plus {1:1.5f} or minus {2:1.5f}".format(*self.L_mcmc))
        print ("We determined C to be {0:1.5f} plus {1:1.5f} or minus {2:1.5f}".format(*self.C_mcmc))

        if self.generateplots:
            print ("Plotting the results in projected space")
            # Plot the results in projected space
            try:
                fig,ax = plt.subplots(figsize=(15,10))
                for opt in samples[np.random.randint(len(samples), size=200)]:
                    ax.plot(self.frequency[self.thresh],\
                            transmissionModel(self.frequency[self.thresh],*opt), color="k", alpha=0.1)
                ax.errorbar(self.frequency[self.thresh][::4], self.ratio_avg[::4],\
                            yerr=self.ratio_err[::4], fmt="b.-",ecolor='b', label='data')
                ax.plot(self.frequency[self.thresh],\
                    transmissionModel(self.frequency[self.thresh], *self.params),\
                    'r-', label='MLE Fit');
                ax.axis('tight');
                ax.legend(loc='best')
                ax.set_xlabel('Frequency [GHz]')
                ax.set_xticklabels(["{0:3.1f}".format(t) for t in ax.get_xticks()])
                ax.set_yticklabels(["{0:1.4f}".format(t) for t in ax.get_yticks()])
                fig.savefig(self.plt_savename + 'projected-space.png')
                plt.close()
            except RuntimeError as err:
                print ("RuntimeError: {0}.\nPlot wasn't generated. Continuing".format(err))
            print ("All done with finding the errorbars! Saving the data now")

    def savedata(self):
        filename = 'mesh_filters.hdf5'
        datafile = h5py.File(filename, 'a')
        try:
            curr_grp = datafile[self.hdf5_name]
        except KeyError:
            curr_grp = datafile.create_group(self.hdf5_name)

        print ("\nSaving data in file {0} under group {1}".format(filename, self.hdf5_name))
        
        savedataset(curr_grp, self.A_mcmc, 'A_mcmc')
        savedataset(curr_grp, self.R_mcmc, 'R_mcmc')
        savedataset(curr_grp, self.T_mcmc, 'T_mcmc')
        savedataset(curr_grp, self.L_mcmc, 'L_mcmc')
        savedataset(curr_grp, self.C_mcmc, 'C_mcmc')
        savedataset(curr_grp, self.etalonlength, 'etalonlength')
        savedataset(curr_grp, self.params, 'params')
        savedataset(curr_grp, self.frequency, 'frequency')
        savedataset(curr_grp, self.thresh, 'thresh')
        savedataset(curr_grp, self.ratio_avg, 'ratio')
        savedataset(curr_grp, self.ratio_err, 'ratio-err')

        try:
            samp_grp = curr_grp['sample']
        except KeyError:
            samp_grp = curr_grp.create_group('sample')

        try:
            nosamp_grp = curr_grp['no-sample']
        except KeyError:
            nosamp_grp = curr_grp.create_group('no-sample')

        try:
            one_grp = curr_grp['onearm']
        except KeyError:
            one_grp = curr_grp.create_group('onearm')

        recursivesave(samp_grp, self.sampledata)
        recursivesave(nosamp_grp, self.nosampledata)
        recursivesave(one_grp, self.onearmdata)

        datafile.close()
        return
