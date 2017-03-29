import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

##
import sys
import random

## equation class

coeffs = [  [4, 5,  6,  7   ],
            [0, 9,  10, 11  ],
            [0, 0,  12, 13  ],
            [0, 0,  0,  14  ]]


def van_der_pol_oscillator_deriv(x, t):
    nx0 = x[1]
    nx1 = -mu * (x[0] ** 2.0 - 1.0) * x[1] - x[0]
    res = np.array([nx0, nx1])
    return res

# ts = np.linspace(0.0, 50.0, 500)

# xs = odeint(van_der_pol_oscillator_deriv, [0.2, 0.2], ts)
# plt.plot(xs[:,0], xs[:,1])
# xs = odeint(van_der_pol_oscillator_deriv, [-3.0, -3.0], ts)
# plt.plot(xs[:,0], xs[:,1])
# xs = odeint(van_der_pol_oscillator_deriv, [4.0, 4.0], ts)
# plt.plot(xs[:,0], xs[:,1])
# plt.gca().set_aspect('equal')
# plt.savefig('vanderpol_oscillator.png')
# plt.show() 



def quadtheta(t, y):


    w_dot = 0.
    x_dot = 0.
    y_dot = 0.
    z_dot = 0.
    return w_dot, x_dot, y_dot, z_dot

## TODO
# define coefficient object, a[i][j] = ...

# define evaluation

class FourQuad:

    def __init__(self):
        coeffs = None

    def evaluate(self, z):
        """
        Evaluates a complex polynomial at z
        """
        total = _zero()
        coeffs = self.coeffs
        for i in range(len(coeffs)):
            total = total + coeffs[i] * z**i
        return total

    def differentiate(self, z):
        return

    def __iter__(self):
        return 

    def __getitem__(self, key):
        return key



## Next steps: 
    ## Need to define class, coefficient calling
        ## Just use scalars, can switch to intervals easily
        ## How general? 

    ## Visualiation of scaling by one dimension
    ## Autodifferentiation

    ## viz to discard big parts of the space, esp. spirals, etc., 
    ## help decide where not to start

def quad_distance(w, x, y, z):
    return [w[i]**2 + x[i]**2 + y[i]**2 + z[i]**2 for i in range(len(w))]

def quad1(w, x, y, z, s=10, r=28, b=2.667):
    w_dot = x*y - b*z
    x_dot = s*(y - x)
    y_dot = r*x - y - x*z
    z_dot = x*y - b*z
    return w_dot, x_dot, y_dot, z_dot


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

def plot_quad(ws, xs, ys, zs, plot_type = 0, txt = ""):    

    if plot_type == 0:
        print("Plotting Double Plot Quad Viz")
        plt.figure(1)

        plt.subplot(2, 1, 1)
        plt.subplots_adjust(top=0.85)
        plt.plot(xs, ws)
        #plt.yscale('linear')
        plt.title('xy')
        plt.grid(True)
        #plt.gca().set_aspect('equal')

        plt.subplot(2, 1, 2)
        plt.plot(ys, zs)
        #plt.yscale('linear')
        plt.title('wz')
        plt.grid(True)
        #plt.gca().set_aspect('equal')
        plt.suptitle(txt, fontsize=14)

        plt.show()

    if plot_type == 1:
        print("Plotting Overlain Double Plot Quad Viz")
        plt.figure(1)

        plt.plot(xs, ws)

        plt.plot(ys, zs)
        #plt.yscale('linear')
        plt.title('x-w, y-z')
        plt.grid(True)
        #plt.gca().set_aspect('equal')
        plt.suptitle(txt, fontsize=14)

        plt.show()

    elif plot_type == 2:
        print("Plotting Sphere Plot Quad Viz")
        fig = plt.figure()
        ax = fig.gca(projection='3d')

        plt.subplots_adjust(top=0.85)
        plt.suptitle(txt, fontsize=14)

        qdist = quad_distance(ws, xs, ys, zs)

        ws = np.divide(ws, qdist)
        xs = np.divide(xs, qdist)
        ys = np.divide(ys, qdist)
        zs = np.divide(zs, qdist)
        
        ax.plot(xs, ys, zs)
        ax.set_xlabel("X Axis")
        ax.set_ylabel("Y Axis")
        ax.set_zlabel("Z Axis")
        ax.set_title("Nonrigorous Solution")
        

        plt.show()

    else: 
        print("Invalid Plot Type")

def main(argv):

    sim = 'demo3'

    if sim == 'demo1':
        dt = 0.01
        stepCnt = 10000

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
        plot_quad(ws, xs, ys, zs, float(argv[1]))

    elif sim == 'demo2':
        dt = 0.01
        stepCnt = 100000

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
            w_dot, x_dot, y_dot, z_dot = quad2(ws[i], xs[i], ys[i], zs[i])
            ws[i + 1] = ws[i] + (w_dot * dt)
            xs[i + 1] = xs[i] + (x_dot * dt)
            ys[i + 1] = ys[i] + (y_dot * dt)
            zs[i + 1] = zs[i] + (z_dot * dt)
        
        plot_quad(ws, xs, ys, zs, float(argv[1]))


    elif sim == 'demo3':
        """
        Loop through simulations
        """
        for i in range(1):

            lambda_1 = .086 # random.random()
            lambda_2 = .141 # (1 - lambda_1) * random.random()
            lambda_3 = 1 - lambda_1 - lambda_2

            dt = 0.01
            stepCnt = 100000

            # Need one more for the initial values
            ws = np.empty((stepCnt + 1,))
            xs = np.empty((stepCnt + 1,))
            ys = np.empty((stepCnt + 1,))
            zs = np.empty((stepCnt + 1,))

            # Setting initial values
            cval = 1

            ws[0], xs[0], ys[0], zs[0] = (  .123146774905 + cval * random.random(), 
                                            0.137604366079 + cval * random.random(), 
                                            -0.0679854330702 + cval * random.random(), 
                                            -0.16021965379  + cval * random.random())

            # Stepping through "time".
            for i in range(stepCnt):
                # Derivatives of the W, X, Y, Z state
                w_dot, x_dot, y_dot, z_dot = quad2(ws[i], xs[i], ys[i], zs[i],
                                                    lambda_1, lambda_2, lambda_3)

                ws[i + 1] = ws[i] + (w_dot * dt)
                xs[i + 1] = xs[i] + (x_dot * dt)
                ys[i + 1] = ys[i] + (y_dot * dt)
                zs[i + 1] = zs[i] + (z_dot * dt)

            # display initial value
            print("w_0, x_0, y_0, z_0 = " 
                                + str(ws[0]) + ", "
                                + str(xs[0]) + ", "
                                + str(ys[0]) + ", "
                                + str(zs[0]))

            # display parameters
            print("lambda_1, lambda_2, lambda_3 = " 
                                + str(lambda_1) + ", "
                                + str(lambda_2) + ", "
                                + str(lambda_3))

            txt = ("Parameters: lambda_1, lambda_2, lambda_3 = " 
                                + str(round(lambda_1, 3)) + ", "
                                + str(round(lambda_2, 3)) + ", "
                                + str(round(lambda_3, 3)) + '\n'
                + "Initial Point: w_0, x_0, y_0, z_0 = " 
                                + str(round(ws[0], 3)) + ", "
                                + str(round(xs[0], 3)) + ", "
                                + str(round(ys[0], 3)) + ", "
                                + str(round(zs[0], 3)) )

            plot_quad(ws, xs, ys, zs, float(argv[1]), txt = txt)

if __name__=="__main__":
    main(sys.argv)


### #

# TODO
# Plot overlain
# plot many areas
