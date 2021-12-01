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

pi_e = (n*q**2/(eps*me))**0.5
omega_e = -q*B/me


w = np.arange(0, 100*abs(omega_e), 0.01*abs(omega_e)) + 0.1

def dispersion(theta, w):
    Xe = (pi_e/w)**2
    Ye = omega_e/w

    R = 1 - Xe/(1 + Ye)
    L = 1 - Xe/(1 - Ye)

    S = (R + L)*0.5
    D = (R - L)*0.5
    P = 1 - Xe

    A = S*(np.sin(theta))**2 + P*(np.cos(theta))**2
    B = R*L*(np.sin(theta))**2 + P*S*(1+(np.cos(theta))**2)
    C = P*R*L
    F = (B**2 - 4*A*C)**0.5

    n1 = 1 - 2*(A - B + C)/(2*A - B + F)
    n2 = 1 - 2*(A - B + C)/(2*A - B - F)
    for i in range(w.size):
        if n1[i] < 0:
            n1[i] = np.nan
    for i in range(w.size):
        if n2[i] < 0:
            n2[i] = np.nan
    
    return n1,n2

n1_0, n2_0 = dispersion(0, w)
n1_30, n2_30 = dispersion(np.pi/6, w)
n1_60, n2_60 = dispersion(np.pi/3, w)
n1_90, n2_90 = dispersion(np.pi/2, w)

k1_0, k2_0 = w/c*n1_0, w/c*n2_0
k1_30, k2_30 = w/c*n1_30, w/c*n2_30
k1_60, k2_60 = w/c*n1_60, w/c*n2_60
k1_90, k2_90 = w/c*n1_90, w/c*n2_90

plt.figure()
plt.plot(k1_0, w/abs(omega_e), label=r"$k+ \theta=0$")
plt.plot(k2_0, w/abs(omega_e), label=r"$k- \theta=0$")
plt.xscale('log')
plt.xlabel('k [/m]')
plt.ylabel(r'$\omega/\Omega_e$')
plt.legend()
plt.show()

plt.figure()
plt.plot(k1_30, w/abs(omega_e), label=r"$k+ \theta=30$")
plt.plot(k2_30, w/abs(omega_e), label=r"$k- \theta=30$")
plt.xscale('log')
plt.xlabel('k [/m]')
plt.ylabel(r'$\omega/\Omega_e$')
plt.legend()
plt.show()

plt.figure()
plt.plot(k1_60, w/abs(omega_e), label=r"$k+ \theta=60$")
plt.plot(k2_60, w/abs(omega_e), label=r"$k- \theta=60$")
plt.xscale('log')
plt.xlabel('k [/m]')
plt.ylabel(r'$\omega/\Omega_e$')
plt.legend()
plt.show()

wL =(omega_e + (omega_e**2 + 4*pi_e**2)**0.5)/2
wR =(-omega_e + (omega_e**2 + 4*pi_e**2)**0.5)/2
w_uh = (pi_e**2 + omega_e**2)**0.5
plt.figure()
plt.plot(k1_90, w/abs(omega_e), label=r"$k+ \theta =90$")
plt.plot(k2_90, w/abs(omega_e), label=r"$k- \theta =90$")
plt.plot(k2_90, w_uh/abs(omega_e)*np.ones(w.size), label="upper hybrid resonance")
plt.plot(k2_90, wR/abs(omega_e)*np.ones(w.size), label="cut off1")
plt.plot(k2_90, wL/abs(omega_e)*np.ones(w.size), label="cut off2")
plt.xscale('log')
plt.xlabel('k [/m]')
plt.ylabel(r'$\omega/\Omega_e$')
plt.legend()
plt.show()