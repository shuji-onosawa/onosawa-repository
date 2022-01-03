from math import pi
import numpy as np
import matplotlib.pyplot as plt
#event1 parameter
ne = 6e+7
nh = ne
q = 1.6e-19
eps = 8.9e-12
myu = 1.3e-6
me = 9.1e-31
mh = 1.7e-27
rho = mh*nh 
B0 = 9e-6
c = 3e+8

pi_e = (ne*q**2/(eps*me))**0.5
pi_h = (nh*q**2/(eps*mh))**0.5
omega_e = -q*B0/me
omega_h = q*B0/mh

def dispersion(theta, w):
    Xe = (pi_e/w)**2
    Xh = (pi_h/w)**2
    Ye = omega_e/w
    Yh = omega_h/w

    R = 1 - Xe/(1 + Ye) - Xh/(1 + Yh) 
    L = 1 - Xe/(1 - Ye) - Xh/(1 - Yh) 

    S = (R + L)*0.5
    D = (R - L)*0.5
    P = 1 - Xe - Xh

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

theta = np.radians(0)
alpha = 50

pi_s = pi_e
omega_s = abs(omega_e)
w = omega_s*np.arange(0.001, 3, 0.001) 

kL1_0, kL2_0, kR1_0, kR2_0, kl1_0, kl2_0 = dispersion(np.degrees(0), w)
kL1_45, kL2_45, kR1_45, kR2_45, kl1_45, kl2_45 = dispersion(np.degrees(45), w)
kL1_90, kL2_90, kR1_90, kR2_90, kl1_90, kl2_90 = dispersion(np.degrees(90), w)

va = B0/(myu*rho)**0.5
wuh = (omega_e**2 + pi_e**2)**0.5
wlh = ((pi_h**2 + omega_h**2)/(1 + pi_e**2/omega_e**2 ))**0.5

plt.figure()
#plt.rcParams["font.size"] = 18
plt.plot(kL1_0, w/abs(omega_s), label = 'L0', color = 'orange')
plt.plot(kL2_0, w/abs(omega_s), color = 'orange')
plt.plot(kL1_45, w/abs(omega_s), label = 'L45', color = 'gold')
plt.plot(kL2_45, w/abs(omega_s), color = 'gold')
plt.plot(kL1_90, w/abs(omega_s), label = 'L90', color = 'red')
plt.plot(kL2_90, w/abs(omega_s), color = 'red')
plt.plot(kR1_0, w/abs(omega_s), label = 'R0', color = 'blue')
plt.plot(kR2_0, w/abs(omega_s), color = 'blue')
plt.plot(kR1_45, w/abs(omega_s), label = 'R45', color = 'slateblue')
plt.plot(kR2_45, w/abs(omega_s), color = 'slateblue')
plt.plot(kR1_90, w/abs(omega_s), label = 'R90', color = 'purple')
plt.plot(kR2_90, w/abs(omega_s), color = 'purple')
plt.hlines(pi_s/abs(omega_s), 0.0001, 1, linestyles="solid", label='fp')
plt.hlines(1, 0.0001, 1, linestyles="solid", colors='k')
plt.hlines(wlh/abs(omega_s), 0.0001, 1, linestyles="dashed", label='lower hybrid f')
plt.hlines(wuh/abs(omega_s), 0.0001, 1, linestyles="dashdot", label='upper hybrid f')
plt.xscale('log')
plt.xlabel('k [/m]')
plt.ylabel(r'$\omega/\Omega_s$')
plt.legend()
plt.show()

print(pi_h/omega_h)