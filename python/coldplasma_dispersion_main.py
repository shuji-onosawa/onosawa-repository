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

theta = 60

omega_s = omega_h
w = omega_s*np.arange(5, 9, 0.0001) 

kL1, kL2, kR1, kR2, kl1, kl2 = dispersion(np.deg2rad(theta), w)
va = B0/(myu*rho)**0.5
wuh = (omega_e**2 + pi_e**2)**0.5
wlh = ((pi_h**2 + pi_he**2 + pi_o**2) / (1 + (pi_e/omega_e)**2 ))**0.5
#kL1, kL2, kR1, kR2, kl1, kl2 = kL1*va/omega_s, kL2*va/omega_s, kR1*va/omega_s, kR2*va/omega_s, kl1*va/omega_s, kl2*va/omega_s

v = (2*q*25/mo)**0.5
omega_s = 2*np.pi
plt.figure()
plt.rcParams["font.size"] = 14
plt.plot(kL1, w/omega_s, label = 'L', color = 'red')
plt.plot(kL2, w/omega_s, color = 'red')
plt.plot(kR1, w/omega_s, label = 'R', color = 'blue')
plt.plot(kR2, w/omega_s, color = 'blue')
plt.hlines(omega_h/omega_s, 0, 1, linestyles='dashed')
plt.hlines(omega_he/omega_s, 0, 1, linestyles='dashed')
plt.hlines(omega_o/omega_s, 0, 1, linestyles='dashed')
#plt.hlines(-omega_e/omega_s, 0, 1, linestyles='dashed')
#plt.hlines(pi_h/omega_s, 0, 1, linestyles='dashed')
#plt.hlines(pi_he/omega_s, 0, 1, linestyles='dashed')
#plt.hlines(pi_o/omega_s, 0, 1, linestyles='dashed')
#plt.hlines(wuh/omega_s, 0, 1, colors='black', linestyles='dashed')
plt.hlines(wlh/omega_s, 0, 1, colors='black', linestyles='dashed')
plt.xscale('log')
plt.xlabel(r'$k [/m]$')
plt.ylabel(r'$\omega [Hz]$')
plt.legend()
plt.show()

print(pi_h/omega_h)
#pitch angle: H+:10-150, O+:20-120
omega = omega_h
m = mh
v10_100 = (w - omega)/kL1/np.cos(np.radians(100 - theta))
v10_110 = (w - omega)/kL1/np.cos(np.radians(110 - theta))
v10_120 = (w - omega)/kL1/np.cos(np.radians(120 - theta))
#v10_130 = (w - omega)/kL1/np.cos(np.radians(130 - theta))
#v10_140 = (w - omega)/kL1/np.cos(np.radians(140 - theta))
#v10_150 = (w - omega)/kL1/np.cos(np.radians(150 - theta))

v20_100 = (w - omega)/kL2/np.cos(np.radians(100 - theta))
v20_110 = (w - omega)/kL2/np.cos(np.radians(110 - theta))
v20_120 = (w - omega)/kL2/np.cos(np.radians(120 - theta))
#v20_130 = (w - omega)/kL2/np.cos(np.radians(130 - theta))
#v20_140 = (w - omega)/kL2/np.cos(np.radians(140 - theta))
#v20_150 = (w - omega)/kL2/np.cos(np.radians(150 - theta))
        
E10_100 = (m*v10_100**2)/2/q
E10_110 = (m*v10_110**2)/2/q
E10_120 = (m*v10_120**2)/2/q
#E10_130 = (m*v10_130**2)/2/q
#E10_140 = (m*v10_140**2)/2/q
#E10_150 = (m*v10_150**2)/2/q

E20_100 = (m*v20_100**2)/2/q
E20_110 = (m*v20_110**2)/2/q
E20_120 = (m*v20_120**2)/2/q
#E20_130 = (m*v20_130**2)/2/q
#E20_140 = (m*v20_140**2)/2/q
#E20_150 = (m*v20_150**2)/2/q

cm = plt.cm.get_cmap('RdYlBu')
plt.figure()
plt.rcParams["font.size"] = 14
plt.plot(E10_100, w/omega_s, label = '100°', color = 'blue')
plt.plot(E20_100, w/omega_s, color = 'blue')
plt.plot(E10_110, w/omega_s, label = '110°', color = 'red')
plt.plot(E20_110, w/omega_s, color = 'red')
plt.plot(E10_120, w/omega_s, label = '120°', color = 'green')
plt.plot(E20_120, w/omega_s, color = 'green')
plt.hlines(omega_h/omega_s, 0, 1e12, colors = 'k', linestyles='solid')
plt.hlines(omega_he/omega_s, 0, 1e12, colors = 'k',linestyles='solid')
plt.hlines(omega_o/omega_s, 0, 1e12, colors = 'k',linestyles='solid')
plt.xlabel('$Energy [eV]$')
plt.ylabel('$\omega [Hz]$')
plt.xscale('log')
plt.legend()
plt.show() 

cm = plt.cm.get_cmap('RdYlBu')


