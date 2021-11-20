import numpy as np
import matplotlib.pyplot as plt

nh = 1e+11
ne = nh 
B = 1e-4
eps = 8.9e-12
me = 9.1e-31
mh = 1.7e-27
mo = 2.7e-26
q = 1.6e-19

pi_e = (ne*q**2/(eps*me))**0.5 / (2*np.pi)
pi_h = (nh*q**2/(eps*mh))**0.5 / (2*np.pi)
omega_e = q*B/me / (2*np.pi)
omega_h = q*B/mh / (2*np.pi) 
omega_o = q*B/mo / (2*np.pi)
wL = (- omega_e + (omega_e**2 + 4*pi_e**2)**0.5)*0.5/pi_e
wR = ( omega_e + (omega_e**2 + 4*pi_e**2)**0.5)*0.5/pi_e

print(pi_e)
print(pi_h)
print(omega_e)
print(omega_h)
print(omega_o)

#w1 = 0.001 + np.arange(1, 3, 0.01)
#k1 = np.power(w1**2 - 1, np.full(w1.size, 0.5))

#w2 = 0.001 + np.arange(wL, 3, 0.01)
#k2 = np.power(w2**2 - w2/(w2 + omega_e/pi_e), np.full(w2.size, 0.5))

#w3 = 0.001 + np.arange(wR, 3, 0.01)
#k3 = np.power(w3**2 - w3/(w3 - omega_e/pi_e), np.full(w3.size, 0.5))


#plt.figure()
#plt.plot(k1, w1)
#plt.plot(k2, w2)
#plt.plot(k3, w3)
plt.xscale('log')
plt.xlim(0.1,4)
plt.ylim(0, 3)
plt.xlabel('$ck/\Pi_e$')
plt.ylabel('$w/\Pi_e$')
plt.show()



