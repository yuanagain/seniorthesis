# Plot of the Lorenz Attractor based on Edward Lorenz's 1963 "Deterministic
# Nonperiodic Flow" publication.
# http://journals.ametsoc.org/doi/abs/10.1175/1520-0469%281963%29020%3C0130%3ADNF%3E2.0.CO%3B2
#
# Note: Because this is a simple non-linear ODE, it would be more easily
#       done using SciPy's ode solver, but this approach depends only
#       upon NumPy.

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def lorenz(x, y, z, s=10, r=13, b=2.667):
    x_dot = s*(y - x)
    y_dot = r*x - y - x*z
    z_dot = x*y - b*z
    return x_dot, y_dot, z_dot


dt = 0.01
stepCnt = 10000

# Need one more for the initial values
xs = np.empty((stepCnt + 1,))
ys = np.empty((stepCnt + 1,))
zs = np.empty((stepCnt + 1,))
ws = np.empty((stepCnt + 1,))

# Setting initial values
xs[0], ys[0], zs[0], ws[0]= (0., 1., 1.05, 0.)

# Stepping through "time".
for i in range(stepCnt):
    # Derivatives of the X, Y, Z state
    x_dot, y_dot, z_dot = lorenz(xs[i], ys[i], zs[i])
    xs[i + 1] = xs[i] + (x_dot * dt)
    ys[i + 1] = ys[i] + (y_dot * dt)
    zs[i + 1] = zs[i] + (z_dot * dt)
    ws[i + 1] = i * dt

#fig = plt.figure()
#ax = fig.gca(projection='3d')

plt.figure(1)

plt.subplot(2, 1, 1)
plt.plot(xs, ys)
plt.yscale('linear')
plt.title('xy')
plt.grid(True)
plt.gca().set_aspect('equal')

plt.subplot(2, 1, 2)
plt.plot(ws, zs)
plt.yscale('linear')
plt.title('wz')
plt.grid(True)
plt.gca().set_aspect('equal')

plt.show()


print(ws[0:10])
print(ys)
print(ws)

#plt.show()
