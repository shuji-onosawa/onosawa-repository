import numpy as np
import matplotlib.pyplot as plt

theta = np.linspace(0 + 10**(-4), 2*np.pi - 10**(-4), 1000)
ganma = 5/3
beta = 0.5

alpha = 0.5*ganma*beta

Vpha = abs(np.cos(theta))
Vphf = np.sqrt(0.5*(1 + alpha + np.sqrt((1 + alpha)**2 - 4*alpha*(np.cos(theta))**2)))
Vphs = np.sqrt(0.5*(1 + alpha - np.sqrt((1 + alpha)**2 - 4*alpha*(np.cos(theta))**2)))

plt.plot(Vpha*np.sin(theta), Vpha*np.cos(theta), 'r', label="Alfven Wave")
plt.plot(Vphf*np.sin(theta), Vphf*np.cos(theta), 'b', label="Fast Wave")
plt.plot(Vphs*np.sin(theta), Vphs*np.cos(theta), 'g', label="Slow Wave")
plt.xlabel(r'$v^{ph}_x/v_A$')
plt.ylabel(r'$v^{ph}_z/v_A$')
plt.legend()

plt.savefig('Freidrichs_diagram_beta%.2f.pdf' % beta)