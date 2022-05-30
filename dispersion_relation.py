import numpy as np
from matplotlib import pyplot as plt

n_e =1e+11
q = 1.6e-19
eps = 8.9e-12
m_e = 9.1e-31
B = 1e-4

pi_e = np.sqrt(n_e*q/(eps*m_e))
omega_e = q*B/m_e

w = np.array(100)
k1 = np.sqrt(w**2 -w*pi_e**2 /(w + omega_e))
k1 = np.sqrt(w**2 -w*pi_e**2 /(w - omega_e))

ax = plt.gca()
plt.figure()
plt.plot(k1, w)
ax.set_xscale('log')  # x軸をlogスケールで描く
plt.show