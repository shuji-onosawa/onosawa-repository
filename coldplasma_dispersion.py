import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import solve

n = 1e+11/25.
q = 1.6e-19
eps = 8.9e-12
me = 9.1e-31
mo = 2.7e-26
mp = 1.7e-27
B = 1e-5
c = 3e+8

pi_e = (n*q**2/(eps*me))**0.5
pi_o = (n*q**2/(eps*mo))**0.5
omega_e = -q*B/me
omega_o = q*B/mo

w = np.arange(1e-5*abs(omega_e), 1e-4*abs(omega_e), 1e-7*abs(omega_e)) + 0.001

def dispersion(theta, w):
    Xe = (pi_e/w)**2
    Xo = (pi_o/w)**2
    Ye = omega_e/w
    Yo = omega_o/w

    R = 1 - Xe/(1 + Ye) - Xo/(1 + Yo)
    L = 1 - Xe/(1 - Ye) - Xo/(1 - Yo)

    S = (R + L)*0.5
    D = (R - L)*0.5
    P = 1 - Xe

    A = S*(np.sin(theta))**2 + P*(np.cos(theta))**2
    B = R*L*(np.sin(theta))**2 + P*S*(1+(np.cos(theta))**2)
    C = P*R*L
    F = np.sqrt(B**2 - 4*A*C)

    n1 = (B + F)/(2*A)
    n2 = (B - F)/(2*A)
    
    for i in range(w.size):
        if n1[i] < 0:
            n1[i] = np.nan
    for i in range(w.size):
        if n2[i] < 0:
            n2[i] = np.nan

    return n1,n2


E_res = mo*v_res**2
n1_0, n2_0 = dispersion(0, w)
k1_0, k2_0 = w/c*np.sqrt(n1_0), w/c*np.sqrt(n2_0)

plt.figure()
plt.plot(k1_0, w/abs(omega_e), label=r"$k+ \theta=0$")
plt.plot(k2_0, w/abs(omega_e), label=r"$k- \theta=0$")
plt.xscale('log')
plt.xlabel('k [/m]')
plt.ylabel(r'$\omega/\Omega_e$')
plt.legend()
plt.show()

