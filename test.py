import numpy as np
import matplotlib.pyplot as plt

a = np.array([1, 2, 3, 4])
b = str(a[0])
print(b)
print('$\Theta$')

alpha1 = np.array([[0, 1, 2, 3],
                   [0, 2, 4, 6]])
for i in range(2):
    alpha1[i,:] = a - i
print(alpha1)

plt.figure()
plt.plot(np.arange(0,4), alpha1[0,:])
plt.plot(np.arange(0,4), alpha1[1,:])
plt.show()

print(alpha1)