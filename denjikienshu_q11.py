import numpy as np

myu = 1.26*10**(-6)
Bz = 11*10**(-9)
R = 71492*10**3
r = 30*R
rho = 2.1*1.66*10**(-21)
T = 9.96*3600
omega = 2*np.pi/T

theta = np.arctan(Bz**2 /(myu*rho*r*r*omega*omega))
rho_t = Bz**2 / (myu*r*r*omega)

print(theta*180/np.pi)
print(rho_t)