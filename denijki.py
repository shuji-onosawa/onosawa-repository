import matplotlib.pyplot as plt
import numpy as np
kp = 5
alpha = 8.9E+4
beta = 45/((1 - 0.159*kp + 0.0093*kp**2)**3)
A = 3*((alpha/2)**(2/3)) * (beta**(1/3))
print(alpha)
print(A)
delta = 0.025
xrange = np.arange(-10, 10, delta)
yrange = np.arange(-10, 10, delta)
X, Y = np.meshgrid(xrange,yrange)
#軸の設定
plt.axis([-10, 10, -10, 10])
plt.gca().set_aspect('equal', adjustable='box')
#描画_
Z = beta*(X**2+Y**2)*Y - A*np.sqrt(X**2+Y**2) + alpha
plt.contour(X, Y, Z, [0])
plt.show()