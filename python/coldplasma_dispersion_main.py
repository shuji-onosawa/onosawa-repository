from math import pi
import numpy as np
import matplotlib.pyplot as plt

ne = 6e+7
nh = 0.46*ne
nhe = 0.11*ne
no = 0.43*ne
q = 1.6e-19
eps = 8.9e-12
me = 9.1e-31
mh = 1.7e-27
mhe = 6.7e-27
mo = 2.7e-26
B0 = 9e-6
c = 3e+8

pi_e = (ne*q**2/(eps*me))**0.5
pi_h = (nh*q**2/(eps*mh))**0.5
pi_he = (nhe*q**2/(eps*mhe))**0.5
pi_o = (no*q**2/(eps*mo))**0.5
omega_e = -q*B0/me
omega_h = q*B0/mh
omega_he = q*B0/mhe
omega_o = q*B0/mo

theta = 0
alpha = 80
w = abs(omega_o)*np.arange(0.001, 20, 0.001) + 0.1

def dispersion(theta, w):
    Xe = (pi_e/w)**2
    Xh = (pi_h/w)**2
    Xhe = (pi_he/w)**2
    Xo = (pi_o/w)**2
    Ye = omega_e/w
    Yh = omega_h/w
    Yhe = omega_he/w
    Yo = omega_o/w

    R = 1 - Xe/(1 + Ye) - Xh/(1 + Yh) - Xhe/(1 + Yhe) - Xo/(1 + Yo)
    L = 1 - Xe/(1 - Ye) - Xh/(1 - Yh) - Xhe/(1 - Yhe) - Xo/(1 - Yo)

    S = (R + L)*0.5
    D = (R - L)*0.5
    P = 1 - Xe - Xh - Xhe - Xo 

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
    
    L1 = np.nan*np.zeros(w.size)
    L2 = np.nan*np.zeros(w.size)
    R1 = np.nan*np.zeros(w.size)
    R2 = np.nan*np.zeros(w.size)
    l1 = np.nan*np.zeros(w.size)
    l2 = np.nan*np.zeros(w.size)

    j = 1j
    right = [0, 0, 0]
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

    return kL1, kL2, kR1, kR2, kl1, kl2

kL1, kL2, kR1, kR2, kl1, kl2 = dispersion(theta, w)

#v = (w + omega_o)/k2/np.cos(np.radians(alpha - theta))
#Eo = mo*v**2/2/q
#cm = plt.cm.get_cmap('RdYlBu')

plt.figure()
plt.plot(kL1, w/abs(omega_o), label = 'L', color = 'orange')
plt.plot(kL2, w/abs(omega_o), color = 'orange')
plt.plot(kR1, w/abs(omega_o), label = 'R', color = 'blue')
plt.plot(kR2, w/abs(omega_o), color = 'blue')
plt.plot(kL1, w/abs(omega_o), label = 'l', color = 'k')
plt.plot(kL1, w/abs(omega_o), color = 'k')
plt.hlines(omega_h/omega_o, 0, 0.002, linestyles='dashed')
plt.hlines(omega_he/omega_o, 0, 0.002, linestyles='dashed')
plt.hlines(omega_o/omega_o, 0, 0.002, linestyles='dashed')
plt.xscale('log')
plt.xlabel('k [/m]')
plt.ylabel(r'$\omega/\Omega_o$')
#plt.colorbar(label='Energy [eV]')
plt.legend()
plt.rcParams["font.size"] = 18
plt.show()

w = abs(omega_e)*np.arange(1e-2, 4, 1e-2) + 0.1
kL1, kL2, kR1, kR2, kl1, kl2 = dispersion(theta, w)

#v = (w + omega_o)/k2/np.cos(np.radians(alpha - theta))
#Eo = mo*v**2/2
#cm = plt.cm.get_cmap('RdYlBu')

Va = B0*c*eps**0.5/(mo*n)**0.5
wuh = (omega_e**2 + pi_e**2)**0.5
wlh = (omega_e**2 *omega_o**2 *(1 + pi_o**2/omega_o**2)/(omega_e**2 + pi_e**2))**0.5
plt.figure()
plt.plot(kL1, w/abs(omega_e), label = 'L', color = 'orange')
plt.plot(kL2, w/abs(omega_e), color = 'orange')
plt.plot(kR1, w/abs(omega_e), label = 'R', color = 'blue')
plt.plot(kR2, w/abs(omega_e), color = 'blue')
plt.plot(kL1, w/abs(omega_e), label = 'l', color = 'k')
plt.plot(kL1, w/abs(omega_e), color = 'k')
plt.hlines(pi_e/abs(omega_e), 0, 1, colors='black', linestyles='dashed')
plt.hlines(1, 0, 1, colors='black', linestyles='dashed')
plt.hlines(wuh/abs(omega_e), 0, 1, colors='black', linestyles='dashed')
plt.xscale('log')
plt.xlabel('k [/m]')
plt.ylabel(r'$\omega/\Omega_e$')
plt.legend()
plt.show()

