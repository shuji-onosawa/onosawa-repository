import numpy as np
import matplotlib.pyplot as plt

n = 5e+10
B0 = 3e-5
eps = 8.9e-12
me = 9.1e-31
mo = 2.7e-26
q = 1.6e-19
c = 3e+8

<<<<<<< HEAD
pi_e = (n*q**2/(eps*me))**0.5
pi_o = (n*q**2/(eps*mo))**0.5
omega_e = -q*B0/me
omega_o = q*B0/mo
=======
w = np.arange(0, 10000*pi_o, 0.01*pi_o)
>>>>>>> 18b93b6054dd6db28cc30c710ee1acd0a58e2f6d

w = abs(omega_o)*np.arange(0.01, 6, 0.01) + 0.1

Xe = (pi_e/w)**2
Xo = (pi_o/w)**2
Ye = omega_e/w
Yo = omega_o/w

R = 1 - Xe/(1 + Ye) - Xo/(1 + Yo)
L = 1 - Xe/(1 - Ye) - Xo/(1 - Yo)

n1 = R
n2 = L

for i in range(w.size):
    if n1[i] < 0:
        n1[i] = np.nan

for i in range(w.size):
    if n2[i] < 0:
        n2[i] = np.nan

k1 = (w*R**0.5)/c
k2 = (w*L**0.5)/c

plt.figure()
plt.plot(k1/c, w/abs(omega_o), label=r'$ n^2 = R $')
plt.plot(k2/c, w/abs(omega_o), label=r'$ n^2 = L $')
plt.xscale('log')
plt.xlabel(r'$k[/m]$')
plt.ylabel(r'$w/\Omega$')
plt.legend()
plt.show()