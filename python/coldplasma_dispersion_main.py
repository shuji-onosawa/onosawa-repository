import numpy as np
import matplotlib.pyplot as plt
from numpy.core.function_base import linspace
#event1 parameter
ne = 6e+7
nh = 0.46*ne
nhe = 0.11*ne
no = 0.43*ne
q = 1.6e-19
eps = 8.9e-12
myu = 1.3e-6

me = 9.1e-31
mh = 1.7e-27
mhe = 6.7e-27
mo = 2.7e-26
rho = mo*no + mh*nh + mhe*nhe
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
    for i in range(w.size):
        s = S[i]
        d = D[i]
        p = P[i]
        N1 = n1[i]
        N2 = n2[i]
        Ex_to_Ey1 = j*d*(p - N1*(np.sin(theta))**2)/(s*p - s*N1*(np.sin(theta))**2 - p*N1*(np.cos(theta))**2)
        Ex_to_Ey2 = j*d*(p - N2*(np.sin(theta))**2)/(s*p - s*N2*(np.sin(theta))**2 - p*N2*(np.cos(theta))**2)
        if (np.angle(Ex_to_Ey1) > 0) and (np.angle(Ex_to_Ey1) < np.pi):
            L1[i] = n1[i]
        if np.angle(Ex_to_Ey1) < 0:
            R1[i] = n1[i]
        if (np.angle(Ex_to_Ey1) == 0) or (np.angle(Ex_to_Ey1) == np.pi):
            l1[i] = n1[i]    
        if (np.angle(Ex_to_Ey2) > 0) and (np.angle(Ex_to_Ey2) < np.pi):
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

def k_energy(w, E, theta, alpha):
    k_energy = (w - omega_o)/((2*E/mo)*np.cos(np.deg2rad(theta))*np.cos(np.deg2rad(alpha)))
    return k_energy

#dispersion calc
omega_s = omega_o
w = omega_s*np.arange(0.1, 1.2, 0.0001) 

va = 1e7
wuh = (omega_e**2 + pi_e**2)**0.5
wlh = ((pi_h**2 + pi_he**2 + pi_o**2) / (1 + (pi_e/omega_e)**2 ))**0.5
kL1_0, kL2_0, kR1_0, kR2_0, kl1, kl2 = dispersion(np.deg2rad(0), w)
#kL1_30, kL2_30, kR1_30, kR2_30, kl1, kl2 = dispersion(np.deg2rad(30), w)
kL1_60, kL2_60, kR1_60, kR2_60, kl1, kl2 = dispersion(np.deg2rad(60), w)

omega_s = 2*np.pi
plt.figure()
plt.rcParams["font.size"] = 14
plt.plot(kL1_0/omega_s, w/omega_s, label = 'L0°', color = 'red')
plt.plot(kL2_0/omega_s, w/omega_s, color = 'red')
#plt.plot(kL1_30/omega_s, w/omega_s, label = 'L30°', color = 'orangered')
#plt.plot(kL2_30/omega_s, w/omega_s, color = 'orangered')
plt.plot(kL1_60/omega_s, w/omega_s, label = 'L60°', color = 'darkorange')
plt.plot(kL2_60/omega_s, w/omega_s, color = 'darkorange')
#plt.plot(kR1_0/omega_s, w/omega_s, label = 'R0°', color = 'blue')
#plt.plot(kR2_0/omega_s, w/omega_s, color = 'blue')
#plt.plot(kR1_30/omega_s, w/omega_s, label = 'R30°', color = 'royalblue')
#plt.plot(kR2_30/omega_s, w/omega_s, color = 'royalblue')
#plt.plot(kR1_60/omega_s, w/omega_s, label = 'R60°', color = 'dodgerblue')
#plt.plot(kR2_60/omega_s, w/omega_s, color = 'dodgerblue')
#plt.plot(w/c/omega_s, w/omega_s, label = 'w = ck/100', linestyle = 'dashed',color = 'gold')
plt.plot(k_energy(w, 1.6e-19, 0, 100)/omega_s, w/omega_s, label = '1 eV, θ=0°', color = 'c')
plt.plot(k_energy(w, 1.6e-19, 60, 100)/omega_s, w/omega_s, label = '1 eV, θ=60°', color = 'm')
#plt.hlines(omega_h/omega_s, 0, 1e-3, linestyles='dashed', colors = 'k')
#plt.hlines(omega_he/omega_s, 0, 1e-3, linestyles='dashed', colors = 'k')
plt.hlines(omega_o/omega_s, 0, 1e-3, linestyles='dashed', colors = 'k')
#plt.hlines(pi_h/omega_s, 0, 1e-3, linestyles='dashed')
#plt.hlines(-omega_e/omega_s, 0, 1, linestyles='dashed')
#plt.hlines(pi_he/omega_s, 0, 1, linestyles='dashed')
#plt.hlines(pi_o/omega_s, 0, 1, linestyles='dashed')
#plt.hlines(wlh/omega_s, 0, 1, colors='black', linestyles='dashed')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$k [/m]$')
plt.ylabel(r'$\omega [Hz]$')
plt.xlim(1e-7, 1e-5)
plt.ylim(3, 11)
plt.legend()
plt.show()

#energy calc
#pitch angle: H+:10-150, O+:20-120
alphaH = np.arange(10, 160, 20)
alphaO = np.arange(20, 130, 20)
v1 = np.nan*np.empty((alphaH.size, w.size))
v2 = v1
E1 = v1
E2 = v1
omega = omega_o
m = mo
alpha = alphaO
for i in range(alpha.size):
    v1[i,:] = (w - omega)/kL1/np.cos(np.radians(alpha[i] - theta))       
    v2[i,:] = (w - omega)/kL2/np.cos(np.radians(alpha[i] - theta)) 
    for j in range(w.size):
        if v1[i,j] > 0:
            E1[i,j] = (m*v1[i,j]**2)/2/q
        if v2[i,j] > 0:
            E2[i,j] = (m*v2[i,j]**2)/2/q

plt.figure()
plt.rcParams["font.size"] = 14
for i in range(alpha.size):
    plt.plot(E1[i,:], w/omega_s, label = str(alpha[i])+'°', color = (1.0 - 0.08*i, 0, 0.08*i))
plt.hlines(omega_h/omega_s, 0, 1e12, colors = 'k', linestyles='dashed')
plt.hlines(omega_he/omega_s, 0, 1e12, colors = 'k',linestyles='dashed')
plt.hlines(omega_o/omega_s, 0, 1e12, colors = 'k',linestyles='dashed')
plt.xlabel('$Energy [eV]$')
plt.ylabel('$\omega [Hz]$')
plt.xscale('log')
plt.yscale('log')
plt.xlim(1, 1e2)
plt.title('$\theta$=' + str(theta) + '°')
plt.legend()
plt.show() 

cm = plt.cm.get_cmap('RdYlBu')


