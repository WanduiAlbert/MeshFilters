import numpy as np
import matplotlib.pyplot as plt
from astropy.constants import c
from scipy.fftpack import rfft
c = c.value*1e-6 #Speed of light in mm/ns
import peakdetect
import scipy.optimize as op
from sys import argv

# Load in a list of the coefficients of the fit. Computed using a matlab fitting routine
coeffs = list(np.loadtxt('fit_coefficients.txt'))
coeffs[-1] = 2*np.pi/coeffs[-1] # I fit for the period. Matlab fits for the frequency. 
assert len(argv)==2, "Enter the name of the file containing the data to be fitted.\n"
try:
    nu, P, Perr = np.loadtxt(argv[1], unpack=True)
except IOError:
    print('Cannot open file %s.'%argv[1])

def transmissionModel(nu, *theta):
    # Fits for the first n terms of a fourier series for a period equal to 1.
    """
    Returns the n term fourier series expansion of a function which has a period of 1 in x and has the fourier
    coefficients listed in theta. The last component of theta is a the fundamental period of the oscillations.
    """
    assert (len(theta)%2)==0, "The length of theta must be even"
    n = len(theta)/2 - 1# Remember n isn't the length of the array.It is the number of terms in the fourier expansion.
    fsr = theta[-1]
    c_k = theta[1:n+1]
    s_k = theta[n+1:-1]
    x = nu/fsr
    index = np.arange(1,n+1)
    fourier_cterm = lambda ck, k: ck*np.cos(2*np.pi*x*k)
    fourier_sterm = lambda sk, k: sk*np.sin(2*np.pi*x*k)
    
    return theta[0] + np.sum(map(fourier_cterm, c_k, index),axis=0)+\
    np.sum(map(fourier_sterm, s_k, index),axis=0)
    
def lnlike(theta, x, y, yerr):
    model = transmissionModel(x, *theta)
    inv_sigma2 = 1.0/(yerr**2)
    return -0.5*(np.sum((y-model)**2*inv_sigma2))

nll = lambda *args: -lnlike(*args) # negative log-likelihood

result = op.minimize(nll, coeffs, args=(nu, P, Perr), method='Powell')
chisq =  lnlike(result['x'],nu,P, Perr )
print result['message'] + " Parameters extracted with a chi-squared value of %4.5f" %chisq

n = len(coeffs)/2 - 1
b0_MLE = result['x'][0]
bk_MLE = (result['x'][1:n+1]**2 + result['x'][n+1:-1]**2)**0.5
freq_MLE = result['x'][-1]

# Plot the fit results vs the actual data to see how good the fit is
# fig, ax = plt.subplots(figsize=(15,10))
# ax.errorbar(nu, P, yerr=Perr, fmt='b.-', ecolor='b');
# ax.plot(nu, transmissionModel(nu, *result["x"]), 'r-');
# ax.grid(which='major')
# ax.axis('tight')
# ax.set_xlabel('Frequency [GHz]')

# Plot the fourier coefficients
index = np.arange(1,n+1)
fig2, ax2 = plt.subplots(figsize=(15,10))
ax2.semilogy(index, bk_MLE/(2*b0_MLE), 'rd-', label='Fourier Fit')
ax2.set_xlabel('k')
ax2.legend(loc='best')
ax2.set_ylabel(r'$\textrm{b}_{\textrm{k}}$');

plt.show()

