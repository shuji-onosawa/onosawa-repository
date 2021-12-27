from math import pi
import numpy as np
import matplotlib.pyplot as plt
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

theta = np.radians(30)
alpha = 50

omega_s = abs(omega_h)
w = omega_s*np.arange(0.001, 10, 0.001) 

kL1, kL2, kR1, kR2, kl1, kl2 = dispersion(theta, w)
va = B0/(myu*rho)**0.5
wuh = (omega_e**2 + pi_e**2)**0.5
wlh = (pi_o**2 / (1 + pi_e**2/omega_e**2 ))**0.5
kL1, kL2, kR1, kR2, kl1, kl2 = kL1*va/omega_s, kL2*va/omega_s, kR1*va/omega_s, kR2*va/omega_s, kl1*va/omega_s, kl2*va/omega_s

plt.figure()
plt.rcParams["font.size"] = 14
plt.plot(kL1, w/omega_s, label = 'L', color = 'orange')
plt.plot(kL2, w/omega_s, color = 'orange')
plt.plot(kR1, w/omega_s, label = 'R', color = 'blue')
plt.plot(kR2, w/omega_s, color = 'blue')
#plt.plot(w/omega_s, w/omega_s, label = 'Alfven wave')
#plt.plot(kL1, w/abs(omega_o), label = 'l', color = 'k')
#plt.plot(kL1, w/abs(omega_o), color = 'k')
#plt.hlines(omega_h/omega_s, 0, 1, linestyles='dashed')
#plt.hlines(omega_he/omega_s, 0, 1, linestyles='dashed')
#plt.hlines(omega_o/omega_s, 0, 1, linestyles='dashed')
#plt.hlines(-omega_e/omega_s, 0, 1, linestyles='dashed')
#plt.hlines(pi_h/omega_s, 0, 1, linestyles='dashed')
#plt.hlines(pi_he/omega_s, 0, 1, linestyles='dashed')
#plt.hlines(pi_o/omega_s, 0, 1, linestyles='dashed')
#plt.hlines(wuh/omega_s, 0, 1, colors='black', linestyles='dashed')
#plt.hlines(wlh/omega_s, 0, 1, colors='black', linestyles='dashed')
plt.xscale('log')
plt.xlabel(r'$kV_a / \Omega_p$')
plt.ylabel(r'$\omega/\Omega_p$')
plt.legend()
plt.show()
""" 
v1 = (w + omega_o)/kL1/np.cos(np.radians(alpha - theta))
v2 = (w + omega_o)/kL2/np.cos(np.radians(alpha - theta))
Eo1 = mo*v1**2/2/q
Eo2 = mo*v2**2/2/q
cm = plt.cm.get_cmap('RdYlBu')
plt.figure()
plt.scatter(kL1, w/omega_o, c=Eo1, cmap=cm)
plt.scatter(kL2, w/omega_o, c=Eo2, cmap=cm)
plt.xlabel('$k [/m]$')
plt.ylabel('$\omega / \Omega_e$')
plt.xscale('log')
plt.colorbar(label='Energy [eV]')
plt.rcParams["font.size"] = 18
plt.show() 

w = 0.999*omega_o*np.ones(1)
kL1, kL2, kR1, kR2, kl1, kl2 = dispersion(theta, w)
v1 = (w + omega_o)/kL1/np.cos(np.radians(alpha - theta))
v2 = (w + omega_o)/kL2/np.cos(np.radians(alpha - theta))
Eo1 = mo*v1**2/2/q
Eo2 = mo*v2**2/2/q
print(Eo1)
print(Eo2) """
