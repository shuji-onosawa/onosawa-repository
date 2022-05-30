import numpy as np
import matplotlib.pyplot as plt


n = 1e+11
q = 1.6e-19
eps = 8.9e-12
me = 9.1e-31
mo = 2.7e-26
mp = 1.7e-27
B = 1e-5
c = 3e+8

pi_e = (n*(-q)**2/(eps*me))**0.5
pi_o = (n*q**2/(eps*mo))**0.5
pi_p = (n*q**2/(eps*mp))**0.5
omega_e = -q*B/me
omega_o = q*B/mo
omega_p = q*B/mp

w = np.arange(0, 10*omega_o, 0.01*omega_o) + 0.1

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

n1_0 = 1 - 2*(A - B + C)/(2*A - B + F)
n2_0 = 1 - 2*(A - B + C)/(2*A - B - F)

for i in range(w.size):
    if n1_0[i] < 0:
        n1_0[i] = np.nan

for i in range(w.size):
    if n2_0[i] < 0:
        n2_0[i] = np.nan

k1_0 = w/c*n1_0**0.5
k2_0 = w/c*n2_0**0.5



plt.figure()
plt.plot(k1_0, w/omega_o, label="k+")
plt.plot(k2_0, w/omega_o, label="k-")
plt.plot(k2_0, np.ones(k2_0.size), label="$\Omega_o$")
plt.xscale('log')
plt.xlabel('k [/m]')
plt.ylabel(r'$\frac{\omega}{\Omega_o}$')
plt.ylim(0,3)
plt.legend()
plt.show()

