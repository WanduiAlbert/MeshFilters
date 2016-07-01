import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from astropy.constants import c
from scipy.fftpack import fft, fftfreq
c = c.value*1e-6 #Speed of light in mm/ns

datadir = '../data/20160412/'
files = glob(datadir + 'fts_150GHz_20160412_0*.txt') 

def quadratic(x, A, B, C):
    return A*x**2 + B*x + C

spectra = []
phase = []
freqs = []
for f in files:
    index, time, encoder, signal = np.loadtxt(f, comments='#', skiprows=17, unpack=True)
    
    # Find the right position of the origin of the interferogram
    peak = np.where(signal == max(signal))[0][0]
    y = signal[peak-1: peak+2]
    x = encoder[peak-1:peak+2]
    X = np.vstack([x**2, x, np.ones_like(x)]).T
    A = np.linalg.solve(X, y)
    actual_peak = -A[1]/(2*A[0])
    sig_peak = quadratic(actual_peak, *list(A))
    encoder -= actual_peak
    encoder *=2 # Encoder now represents the path length difference
    
#     peak = encoder[np.where(signal == max(signal))[0]][0]
#     thresh = np.logical_and(encoder > (peak - 1.0), encoder < (peak + 1.0))
#     encoder = encoder[thresh]
#     signal = signal[thresh]

    # Remove the dc component as much as we can
    mean = np.average(signal)
    signal -= mean
    
    y = fft(signal)
    d = np.mean(np.diff(encoder))
    n = len(encoder)
    k = fftfreq(n,d)*c
    spectra += [np.abs(y)]
    phase += [np.angle(y)]
    freqs += [k]
spectra = np.array(spectra)
phase = np.array(phase)
freqs = np.array(freqs)

ypol = np.arange(7)
xpol = np.arange(7, 14)

# Average out our spectra to reduce the noise
pol = xpol
code = 'x'
nosample = np.average(spectra[pol][::2], axis=0)
nosamplephase = np.average(phase[pol][::2], axis=0)
nosampleerr = np.std(spectra[pol][::2], axis=0)

sample = np.average(spectra[pol][1::2], axis=0)
samplephase = np.average(phase[pol][::2], axis=0)
sampleerr = np.std(spectra[pol][1::2], axis=0)

k = np.average(freqs[pol], axis=0)
kerr = np.std(freqs[pol], axis=0)

thresh = np.logical_and(k > 138, k < 160)

fig, ax = plt.subplots(figsize=(15,10))
ax.plot(k[thresh], nosample[thresh],'k', label='Reference Spectra');
ax.plot(k[thresh], sample[thresh],'r', label='Sample Spectra');
ax.axis('tight');
ax.legend(loc='best')
ax.set_xlabel('Frequency, GHz');
plt.savefig('spectra.png')
plt.show()

transmission = sample/nosample
t_err = transmission * np.sqrt((sampleerr/sample)**2 + (nosampleerr/nosample)**2)
fig, ax =plt.subplots(figsize=(15,10))
ax.errorbar(k[thresh], transmission[thresh], yerr = t_err[thresh], fmt = 'b-o', ecolor='r')
ax.grid(which='major')
ax.set_xlabel('Frequency, GHz');
ax.set_ylabel('ratio')
ax.axis('tight');
ax.set_ylim([0.8, 1.2])
plt.savefig('transmission.png')
plt.show()

#Let's try some interpolation
from scipy import interpolate
f = interpolate.interp1d(k[thresh], transmission[thresh], 'cubic')
kmin = k[thresh][0]
kmax = k[thresh][-1]
k_new = np.r_[kmin:kmax:2*len(k)*1j]
t_new = f(k_new)

fig, ax =plt.subplots(figsize=(15,10))
ax.plot(k_new, t_new,'r', label= 'interpolated FP spectrum')
# ax.errorbar(k[thresh], transmission[thresh], yerr = t_err[thresh], color = 'b')
ax.grid(which='major')
ax.set_xlabel('Frequency, GHz');
ax.set_ylabel('Transmission')
ax.axis('tight');
from peakdetect import peakdetect
max_peaks, min_peaks = peakdetect(t_new, k_new)
loc, val = map(list, zip(*max_peaks))
ax.plot(loc, val, 'bd', markersize=5, label='peaks');
ax.legend(loc='best');
plt.savefig('interpolated.png')
plt.show()

L = 2.25*25.4 # thickness of the etalon in mm
from scipy.optimize import curve_fit

def get_transmission(nu, y, F, A, B):
    return y*(1 - F * (1 - np.cos(A*nu + B))/2)

R = 0.00207
F = 4*R/(1-R)**2
A = 4*np.pi*L/c
B = np.pi
y = np.average(val)
guesses = [y, F, A, B]
popt, pcov = curve_fit(get_transmission, k[thresh], transmission[thresh], p0=guesses, sigma=t_err[thresh])
# popt, pcov = curve_fit(get_transmission, k_new, t_new, p0=guesses)
perr =np.sqrt(np.diag(pcov))

fig, ax =plt.subplots(figsize=(15,10))
# ax.plot(k[thresh], transmission[thresh], 'r-', markersize=10)
ax.plot(k_new, t_new,'r', label='interpolated data')
# ax.errorbar(k[thresh], transmission[thresh], yerr = t_err[thresh], fmt = 'r-o',label='data')
ax.plot(k_new, get_transmission(k_new, *popt), 'b', label='fit')
ax.grid(which='major')
ax.set_xlabel('Frequency, GHz');
ax.set_ylabel('Transmission')
ax.axis('tight');
ax.legend(loc='best');
plt.savefig('interpolatedfit.png')
plt.show()

# Calculate our relevant parameters and errors on them
y, F, a, b = popt[0], popt[1], popt[2], popt[3]
sigy, sigF, siga, sigb = perr[0], perr[1], perr[2], perr[3]
A = (1 - F/4)*(1 - np.sqrt(y))
T = (1 - F/4)*(np.sqrt(y))
R = (2 + F - 2*(1+F)**0.5)/F
sigA = A*np.sqrt((sigy/(2*y**0.5*(1-y**0.5)))**2 + (sigF/(4*(1-F/4)))**2)
sigT = T*np.sqrt((sigy/(2*y))**2 + (sigF/(4*(1-F/4)))**2)
sigR = R/(F*(1+F)**0.5)*sigF

print "absorption: ", A, sigA
print "Transmission: ", T, sigT
print "Reflection: ", R, sigR

plt.show()
