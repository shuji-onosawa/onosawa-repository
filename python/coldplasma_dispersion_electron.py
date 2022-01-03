import numpy as np
import matplotlib.pyplot as plt

<<<<<<< HEAD
n = 1e+11/25.
=======
n = 0.3e+9
>>>>>>> 44c42c3b706c54ec05d8897cae8f39ba70887dc3
q = 1.6e-19
eps = 8.9e-12
me = 9.1e-31
mo = 2.7e-26
mp = 1.7e-27
B = 1e-5
c = 3e+8

pi_e = (n*q**2/(eps*me))**0.5
omega_e = -q*B/me

theta = np.degrees(60)
w = np.arange(0, 3*abs(omega_e), 0.001*abs(omega_e))

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

    n1 = (B + F)/(2*A)
    n2 = (B - F)/(2*A)
    for i in range(w.size):
        if n1[i] < 0:
            n1[i] = np.nan
    for i in range(w.size):
        if n2[i] < 0:
            n2[i] = np.nan
    L1 = np.nan*np.zeros(w.size)
    L2 = np.nan*np.zeros(w.size)
    R1 = np.nan*np.zeros(w.size)
    R2 = np.nan*np.zeros(w.size)
    l1 = np.nan*np.zeros(w.size)
    l2 = np.nan*np.zeros(w.size)

    j = 1j
    for i in range(w.size):
        s = S[i]
        d = D[i]
        p = P[i]
        N1 = n1[i]
        N2 = n2[i]
        Ex_to_Ey1 = j*d*(p - N1*(np.sin(theta))**2)/(s*p - s*N1*(np.sin(theta))**2 - p*N1*(np.cos(theta))**2)
        Ex_to_Ey2 = j*d*(p - N2*(np.sin(theta))**2)/(s*p - s*N2*(np.sin(theta))**2 - p*N2*(np.cos(theta))**2)
        if np.angle(Ex_to_Ey1) > 0:
            L1[i] = n1[i]
        if np.angle(Ex_to_Ey1) < 0:
            R1[i] = n1[i]
        if (np.angle(Ex_to_Ey1) == 0) or (np.angle(Ex_to_Ey1) == np.pi):
            l1[i] = n1[i]    
        if np.angle(Ex_to_Ey2) > 0:
            L2[i] = n2[i]
        if np.angle(Ex_to_Ey2) < 0:
            R2[i] = n2[i]
        if (np.angle(Ex_to_Ey2) == 0) or (np.angle(Ex_to_Ey2) == np.pi):
            l2[i] = n2[i]        

    kL1 = w/c*np.sqrt(L1)
    kL2 = w/c*np.sqrt(L2)
    kR1 = w/c*np.sqrt(R1)
    kR2 = w/c*np.sqrt(R2)
    kl1 = w/c*np.sqrt(l1)
    kl2 = w/c*np.sqrt(l2)

<<<<<<< HEAD
plt.figure()
plt.plot(k1_0, w/abs(omega_e), label=r"$k+ \theta=0$")
plt.plot(k2_0, w/abs(omega_e), label=r"$k- \theta=0$")
plt.plot(k2_90, pi_e/abs(omega_e)*np.ones(w.size), linestyle='--',label="fp")
plt.xscale('log')
plt.xlabel('k [/m]')
plt.ylabel(r'$\omega/\Omega_e$')
plt.legend()
plt.show()
=======
    return kL1, kL2, kR1, kR2, kl1, kl2
    
kL1_0, kL2_0, kR1_0, kR2_0, kl1_0, kl2_0 = dispersion(np.degrees(0), w)
kL1_45, kL2_45, kR1_45, kR2_45, kl1_45, kl2_45 = dispersion(np.degrees(45), w)
kL1_90, kL2_90, kR1_90, kR2_90, kl1_90, kl2_90 = dispersion(np.degrees(90), w)
wL =(omega_e + (omega_e**2 + 4*pi_e**2)**0.5)/2
wR =(-omega_e + (omega_e**2 + 4*pi_e**2)**0.5)/2
w_uh = (pi_e**2 + omega_e**2)**0.5
>>>>>>> 44c42c3b706c54ec05d8897cae8f39ba70887dc3

print(pi_e/omega_e)
print(wL/omega_e)
print(wR/omega_e)
print(w_uh/omega_e)

plt.figure()
#plt.rcParams["font.size"] = 18
plt.plot(kL1_0, w/abs(omega_e), label = 'L0', color = 'orange')
plt.plot(kL2_0, w/abs(omega_e), color = 'orange')
plt.plot(kL1_45, w/abs(omega_e), label = 'L45', color = 'gold')
plt.plot(kL2_45, w/abs(omega_e), color = 'gold')
plt.plot(kL1_90, w/abs(omega_e), label = 'L90', color = 'red')
plt.plot(kL2_90, w/abs(omega_e), color = 'red')
plt.plot(kR1_0, w/abs(omega_e), label = 'R0', color = 'blue')
plt.plot(kR2_0, w/abs(omega_e), color = 'blue')
plt.plot(kR1_45, w/abs(omega_e), label = 'R45', color = 'slateblue')
plt.plot(kR2_45, w/abs(omega_e), color = 'slateblue')
plt.plot(kR1_90, w/abs(omega_e), label = 'R90', color = 'purple')
plt.plot(kR2_90, w/abs(omega_e), color = 'purple')
plt.hlines(pi_e/abs(omega_e), 0.0001, 1, linestyles="solid", label='fp')
plt.hlines(1, 0.0001, 1, linestyles="solid", colors='k')
plt.hlines(wL/abs(omega_e), 0.0001, 1, linestyles="dashed", label='L cutoff')
plt.hlines(wR/abs(omega_e), 0.0001, 1, linestyles="dashdot", label='R cutoff')
plt.hlines(w_uh/abs(omega_e), 0.0001, 1, linestyles="dotted", label='upper hybrid f')
plt.xscale('log')
plt.xlabel('k [/m]')
plt.ylabel(r'$\omega/\Omega_e$')
plt.legend()
plt.show()

<<<<<<< HEAD
wL =(omega_e + (omega_e**2 + 4*pi_e**2)**0.5)/2
wR =(-omega_e + (omega_e**2 + 4*pi_e**2)**0.5)/2
w_uh = (pi_e**2 + omega_e**2)**0.5
plt.figure()
plt.plot(k1_90, w/abs(omega_e), label=r"$k+ \theta =90$")
plt.plot(k2_90, w/abs(omega_e), label=r"$k- \theta =90$")
plt.plot(k2_90, w_uh/abs(omega_e)*np.ones(w.size), linestyle='--',label="upper hybrid resonance")
plt.plot(k2_90, pi_e/abs(omega_e)*np.ones(w.size), linestyle='--',label="fp")
plt.plot(k2_90, wR/abs(omega_e)*np.ones(w.size), linestyle='--',label="cut off1")
plt.plot(k2_90, wL/abs(omega_e)*np.ones(w.size), linestyle='--',label="cut off2")
plt.xscale('log')
plt.xlabel('k [/m]')
plt.ylabel(r'$\omega/\Omega_e$')
plt.legend()
plt.show()
=======


>>>>>>> 44c42c3b706c54ec05d8897cae8f39ba70887dc3
