import numpy as np
import matplotlib.pyplot as plt

n = 5e+10
q = 1.6e-19
eps = 8.9e-12
me = 9.1e-31
mo = 2.7e-26
mp = 1.7e-27
B0 = 3e-5
c = 3e+8

pi_e = (n*q**2/(eps*me))**0.5
pi_o = (n*q**2/(eps*mo))**0.5
omega_e = -q*B0/me
omega_o = q*B0/mo


def dispersion(theta, w):
    Xe = (pi_e/w)**2
    Xo = (pi_o/w)**2
    Ye = omega_e/w
    Yo = omega_o/w

    R = 1 - Xe/(1 + Ye) - Xo/(1 + Yo)
    L = 1 - Xe/(1 - Ye) - Xo/(1 - Yo)

    S = (R + L)*0.5
    D = (R - L)*0.5
    P = 1 - Xe - Xo

    A = S*(np.sin(theta))**2 + P*(np.cos(theta))**2
    B = R*L*(np.sin(theta))**2 + P*S*(1+(np.cos(theta))**2)
    C = P*R*L
    F = (B**2 - 4*A*C)**0.5

    n1 = (B + F)/(2*A)
    n2 = (B - F)/(2*A)
    for i in range(w.size):
        if n1[i] < 0:
            n1[i] = np.nan
    for i in range(w.size):
        if n2[i] < 0:
            n2[i] = np.nan
    
    G = B**2 - 4*A*C
    for i in range(w.size):
        if G[i] < 0:
            n1[i], n2[i] = n2[i], n1[i]
    
    n1 = n1**0.5
    n2 = n2**0.5

    k1, k2 = w/c*n1, w/c*n2
    return k1,k2

theta = 0
w = abs(omega_e)*np.arange(0.01, 10, 0.01) + 0.1
k1, k2 = dispersion(theta, w)

Va = B0*c*eps**0.5/(mo*n)**0.5
wuh = (omega_e**2 + pi_e**2)**0.5
wlh = (omega_e**2 *omega_o**2 *(1 + pi_o**2/omega_o**2)/(omega_e**2 + pi_e**2))**0.5
plt.figure()
plt.plot(k1, w/abs(omega_e), label=r'$k+$')
plt.plot(k2, w/abs(omega_e), label=r'$k-$')
plt.hlines(pi_e/abs(omega_e), 0, 1, colors='black', linestyles='dashed')
plt.hlines(1, 0, 1, colors='black', linestyles='dashed')
plt.hlines(wuh/abs(omega_e), 0, 1, colors='black', linestyles='dashed')
plt.xscale('log')
plt.xlabel('k [/m]')
plt.ylabel(r'$\omega/\Omega_e$')
plt.legend()
plt.show()

w = abs(omega_o)*np.arange(0.01, 30000, 0.01) + 0.1
k1, k2 = dispersion(theta, w)

plt.figure()
plt.plot(k1, w/abs(omega_o), label=r'$k+$')
plt.plot(k2, w/abs(omega_o), label=r'$k-$')
plt.plot(w/Va, w/abs(omega_o), linestyle='dashed')
plt.hlines(omega_e/omega_o, 0, 0.002, colors='black', linestyles='dashed')
#plt.hlines(wlh/abs(omega_o), 0, 1, colors='black', linestyles='dashed')
plt.xscale('log')
plt.xlabel('k [/m]')
plt.ylabel(r'$\omega/\Omega_o$')
plt.legend()
plt.show()