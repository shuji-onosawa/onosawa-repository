import numpy as np
#for文
myu = [1.26e-6, 1,26e-6]
rho = [5*1.66e-21, 0.5*1.66e-21]
v = [5e5, 5e5]
Beq = [4e-5, 4e-4]
Rp = [6371e3, 71492e3]
Rm =[0, 0]

for i in range(2):
 Rm[i] = Rp[i]*np.sqrt(np.sqrt(Beq[i]**2/(myu[i]*rho[i]*v[i]**2)))
 print(Rm[i]*1e-3, 'km')
 print(Rm[i]/Rp[i])

#zip function
myus = myu
rhos = rho
vs = v
Beqs = Beq
Rps = Rp
Rms = [0, 0]

for myu, rho, v, Beq, Rp, Rm in zip(myus, rhos, vs, Beqs, Rps, Rms):
    Rm = Rp*np.sqrt(np.sqrt(Beq**2/(myu*rho*v**2)))
    print(Rm*1e-3, 'km')
    print(Rm/Rp)

