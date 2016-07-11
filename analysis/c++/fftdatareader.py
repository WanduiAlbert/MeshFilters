import numpy as np
import matplotlib.pyplot as plt
from astropy.constants import c
c = c.value*1e-6 #Speed of light in mm/ns

R = 0.0234
c_k, s_k = np.loadtxt('fft_results.dat', unpack=True)

b0 = c_k[0]
c_k = c_k[1:11]
s_k = s_k[1:11]
b_k = (c_k**2 + s_k**2)**0.5

print "I'm here. Starting to run into trouble!!!!!"

index = np.arange(1, 11)
fig, ax = plt.subplots(figsize=(15,10))
ax.semilogy(index, b_k/(2*b0),'rd-', label='Fourier Series')
ax.semilogy(index, R**index,'bo-', label = 'Actual Series')
ax.set_xlabel('k')
ax.legend(loc='best')
ax.set_ylabel(r'$\textrm{b}_{\textrm{k}}$')


# fig, ax = plt.subplots(figsize=(15,10))
# ax.semilogy(index, b_k/(2*b0)-R**index,'rd-')
# ax.set_xlabel('k')
# ax.set_ylabel(r'$\Delta \textrm{b}_{\textrm{k}}$')
plt.show()