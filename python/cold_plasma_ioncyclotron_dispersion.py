import numpy as np
import matplotlib.pyplot as plt

no = 1e+11
ne = no 
B = 1e-5
eps = 8.9e-12
me = 9.1e-31
mo = 2.7e-26
q = 1.6e-19
c = 3e+8
pi_e = (ne*q**2/(eps*me))**0.5 / (2*np.pi)
pi_o = (no*q**2/(eps*mo))**0.5 / (2*np.pi)
omega_e = q*B/me / (2*np.pi)
omega_o = q*B/mo / (2*np.pi)

w = np.arange(0,10,0.1)

k = (1 - pi_e**2 *(1 + omega_o/omega_e)/((w + omega_e)*(w - omega_o)))*w**2

for i in range(100):
    if k[i] < 0:
        k[i] = np.nan 

k2 = np.power(k, np.full(k.size, 0.5))

plt.figure()
plt.plot(k/c, w)
plt.xscale('log')
""" plt.xlim(0.1,4)
plt.ylim(0, 3) """ 
plt.xlabel('$k[/m]$')
plt.ylabel('$w [Hz]$')
plt.show()