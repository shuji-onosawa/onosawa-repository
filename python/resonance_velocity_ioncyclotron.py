import numpy as np
import sympy as sp
import matplotlib.pyplot as plt


n = 1e+11
q = 1.6e-19
eps = 8.9e-12
me = 9.1e-31
mo = 2.7e-26
mp = 1.7e-27
B = 1e-5
c = 3e+8

pi_e = (n*q**2/(eps*me))**0.5
pi_o = (n*q**2/(eps*mo))**0.5
pi_p = (n*q**2/(eps*mp))**0.5
omega_e = q*B/me
omega_o = q*B/mo
omega_p = q*B/mp

<<<<<<< HEAD:python/resonance_velocity_ioncyclotron.py
w = np.arange(0.99*omega_o, omega_o, 0.001*omega_o) 
=======
w = np.arange(0, 10000*omega_e, 0.001*omega_e) + 0.001
>>>>>>> 18b93b6054dd6db28cc30c710ee1acd0a58e2f6d:python/coldplasma_dispersion.py

Xe = (pi_e/w)**2
Xo = (pi_o/w)**2
Xp = (pi_p/w)**2
Ye = omega_e/w
Yo = omega_o/w
Yp = omega_p/w

R = 1 - Xe/(1 + Ye) - Xo/(1 + Yo)
L = 1 - Xe/(1 - Ye) - Xo/(1 - Yo)

S = (R + L)*0.5
D = (R - L)*0.5
P = 1 - Xe - Xo

theta = 0

A = S*(np.sin(theta))**2 + P*(np.cos(theta))**2
B = R*L*(np.sin(theta))**2 + P*S*(1+(np.cos(theta))**2)
C = P*R*L
F = (B**2 - 4*A*C)**0.5

n2 = 1 - 2*(A - B + C)/(2*A - B - F)

k2 = w/c*n2**0.5
alpha = 80
v = (w + omega_o)/k2/np.cos(np.radians(alpha - theta))
Eo = mo*v**2/2
cm = plt.cm.get_cmap('RdYlBu')
plt.figure()
plt.scatter(k2, w/omega_o, c=v, cmap=cm)
plt.xlabel('$k [/m]$')
plt.ylabel('$\omega / \Omega_e$')
plt.xscale('log')
plt.colorbar(label='v [m/s]')
plt.show()

plt.figure()
plt.scatter(k2, w/omega_o, c=Eo, cmap=cm)
plt.xlabel('$k [/m]$')
plt.ylabel('$\omega / \Omega_e$')
plt.xscale('log')
plt.colorbar(label='Energy [eV]')
plt.show()




