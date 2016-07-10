import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from astropy.constants import c
c = c.value*1e-6 #Speed of light in mm/ns

R = 0.0234
c_k, s_k = np.loadtxt('fft_results.dat', unpack=True)

print c_k[0], s_k[0]

b0 = c_k[0]
c_k = c_k[1:-1]
s_k = s_k[1:-1]
b_k = (c_k**2 + s_k**2)

index = np.arange(1, len(b_k)+1)
print len(index) == len(b_k)
fig, ax = plt.subplots(figsize=(15,10))
ax.plot(index[:10], b_k[:10]/(2*b0),'rd-', label='Fourier Series')
ax.plot(index[:10], R**index[:10],'bo-', label = 'Actual Series')
ax.set_xlabel('k')
ax.legend(loc='best')
ax.set_ylabel(r'$\textrm{b}_{\textrm{k}}$')
plt.show()