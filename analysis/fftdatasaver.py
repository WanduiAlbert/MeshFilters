import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from astropy.constants import c
c = c.value*1e-6 #Speed of light in mm/ns

T = 0.9763
R = 0.0234
A = 0.0003
L = 2.25*25.4
nu_fsr = c/(2*L)

min_x =  int(150/nu_fsr) + 1
max_x = min_x + 1
n = 2**20
x = np.linspace(min_x, max_x, n)
nu = x*nu_fsr
# nu = np.linspace(150, 150 + nu_fsr, n) # Range of frequencies of interest
# x = nu/nu_fsr

F = 4*R/(1-R)**2
delta = 2*np.pi*x
P = (T/(1-R))**2*1./(1 + F*np.sin(delta/2)**2)
# P = 2.0*np.sin(2*np.pi*x)
# mask = np.logical_and(x >= x[n/2], x <= (x[n/2]+1)) #Pick a range that covers one entire period of the data

filename = 'fft_data.dat'
f = open(filename, 'w')

# fig, ax = plt.subplots(figsize=(15,10))
# ax.plot(nu, P)
# ax.set_xlabel(r'Frequency [GHz]')
# ax.set_ylabel(r'Transmission')
# plt.show()


f.write('%d\n' %n)
np.savetxt(f, P)

f.close()

