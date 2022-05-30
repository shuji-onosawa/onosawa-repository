import numpy as np

myu = np.array([1.26e-6,1.26e-6])
Beq = np.array([4e-5, 4e-4])
Bsw = np.array([5e-9, 5e-10])
v = np.array([5e+5,5e+5])
Rp = np.array([6371e+3, 71492e+3])
rho = np.array([5, 0.5])*1.66e-21
Tp = np.array([24, 9.93])*3600
omega = 2*np.pi/Tp
eps = np.array([1e-3, 1e-3])

Rm = Rp*(Beq**2/(myu*rho*v**2))**(1/6)
Ls = (Beq*omega*Rp**2/(4*eps*v*Bsw*Rm))**0.33

print('Rm', Rm)
print('Rm/Rp', Rm/Rp)
print('Ls', Ls)
