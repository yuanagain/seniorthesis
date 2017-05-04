from ad import adnumber
from ad.admath import *


class Newton: 
    def __init__(self, start, poly):
        self.start = start
        self.poly = poly
        self.iterates = 0
        self.deriv = poly.derive()
        self.step = start

    def iterate(self):
        """
        Performs one Newton iteration, returns change between values.
        """
        self.iterates += 1
        x = self.step.midpoint()
        fx = self.poly(x)

        ## iterate on derivative
        ## self.deriv = self.deriv.derive()

        self.step = x - (fx / self.deriv(x))

        ## return the change
        diff = x - self.step
        return diff

    def iterate_until(self, res = 0, max_iterates = 20):
        """
        Iterates until at resolution or until maximum number
        of iterations has been reached. Returns True if convergence
        achieved, returns False otherwise.
        """

        while (self.iterates < max_iterates - 1):
            if abs(self.iterate()) < res:
                return True

        if abs(self.iterate()) < res:
                return True

        return False

    def __str__(self):
        """
        Returns string representation
        """
        return "Newton's Iterator\n" + "Start: " + str(self.start) + "\nFunction: " + str(self.poly)



# Goal: find periodicity in solution to ODE

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

## setting up parameters
default_lambda_1 = .2523718
default_lambda_2 = .392931
default_lambda_3 = 1 - default_lambda_1 - default_lambda_2

def quad2(x_1, y_1, x_2, y_2, 
            lambda_1 = default_lambda_1, 
            lambda_2 = default_lambda_2, 
            lambda_3 = default_lambda_3):
    """ 
    dz1/dt = lambda_2 * z1^2 - (lambda_2 + lambda_3) * z1 * z2
    dz2/dt = lambda_1 * z2^2 - (lambda_1 + lambda_3) * z1 * z2
    http://www.math.kit.edu/iag3/~herrlich/seite/wws-11/media/wws-talk-valdez.pdf
    """

    x_1_dot = lambda_2 * (x_1**2 - y_1**2) - (lambda_2 + lambda_3) * (x_1*x_2 - y_1*y_2)
    y_1_dot = 2 * lambda_2 * x_1 * y_1 - (lambda_2 + lambda_3) * (x_1*y_2 + y_1*x_2)
    x_2_dot = lambda_1 * (x_2**2 - y_2**2) - (lambda_1 + lambda_3) * (x_1*x_2 - y_1*y_2)
    y_2_dot = 2 * lambda_1 * x_2 * y_2 - (lambda_1 +lambda_3) * (x_1*y_2 + y_1*x_2)

    return x_1_dot, y_1_dot, x_2_dot, y_2_dot


# For (1), we could try two things. Fix a period T, say 1, then let
# g(x) = f(x,T), where f(x,0) = x. In other words, we integrate a 
# time T from the initial condition x. Then we minimize |g(x)-x| 
# over all x (for example using Newton's method). Another 
# (essentially similar) possibility is to minimize |g(x)-x|^{2}.

## Q: How do we parameterize then?
## Q: Autodiff or minimize tools?

## Write as f(x, T)


import scipy.integrate as integrate


def g(x, T = 1):
    return f(x, T)

def minimize_gdist(x_0):
    """
    Minimize |g(x)-x|
    """

# Setting initial values
xs[0], ys[0], zs[0], ws[0]= (0., 1., 1.05, 0.)

def f(x, T, dt = 0.001):
    stepCnt = round(T / dt)

    # Need one more for the initial values
    ws = np.empty((stepCnt + 1,))
    xs = np.empty((stepCnt + 1,))
    ys = np.empty((stepCnt + 1,))
    zs = np.empty((stepCnt + 1,))

    # Setting initial values
    ws[0], xs[0], ys[0], zs[0] = (  0.372854105052, 
                                    0.393518965248, 
                                    -0.0359026080443, 
                                    -0.216701666067 )
    # Stepping through "time".
    for i in range(stepCnt):
        # Derivatives of the W, X, Y, Z state
        w_dot, x_dot, y_dot, z_dot = quad1(ws[i], xs[i], ys[i], zs[i])
        ws[i + 1] = ws[i] + (w_dot * dt)
        xs[i + 1] = xs[i] + (x_dot * dt)
        ys[i + 1] = ys[i] + (y_dot * dt)
        zs[i + 1] = zs[i] + (z_dot * dt)

    # return the values

    return ws, xs, ys, zs

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
