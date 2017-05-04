import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import numdifftools as nd


default_lambda_1, default_lambda_2, default_lambda_3 = 0.086, 0.141, 0.773
default_start = (0.372854105052, 0.393518965248, -0.0359026080443, -0.216701666067)
x_0 = default_start
res = 0.01
dt = res

def quad_sq_distance(x, y):
    """Computes the squared distance"""
    dists = [ x[i] - y[i] for i in range(len(x) )]
    dists = [ dists[i]**2 for i in range(len(x) )]
    return sum(dists)


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

    elif plot_type == 1:
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

def experiment_1(start_pt = default_start, 
                 T = 1, 
                 lmbda = [default_lambda_1, default_lambda_2, default_lambda_3], 
                 res = 0.001,
                 expmt = "search"):
    
    ## define evaluation function
    def dots(x_0, lmbda):
        """
        dz1/dt = lambda_2 * z1^2 - (lambda_2 + lambda_3) * z1 * z2
        dz2/dt = lambda_1 * z2^2 - (lambda_1 + lambda_3) * z1 * z2
        http://www.math.kit.edu/iag3/~herrlich/seite/wws-11/media/wws-talk-valdez.pdf
        """
        x_1 = x_0[0]
        y_1 = x_0[1]
        x_2 = x_0[2]
        y_2 = x_0[3]

        # print(lmbda)
        lambda_1 = lmbda[0]
        lambda_2 = lmbda[1]
        lambda_3 = lmbda[2]

        x_1_dot = lambda_2 * (x_1**2 - y_1**2) - (lambda_2 + lambda_3) * (x_1*x_2 - y_1*y_2)
        y_1_dot = 2 * lambda_2 * x_1 * y_1 - (lambda_2 + lambda_3) * (x_1*y_2 + y_1*x_2)
        x_2_dot = lambda_1 * (x_2**2 - y_2**2) - (lambda_1 + lambda_3) * (x_1*x_2 - y_1*y_2)
        y_2_dot = 2 * lambda_1 * x_2 * y_2 - (lambda_1 +lambda_3) * (x_1*y_2 + y_1*x_2)

        return [x_1_dot, y_1_dot, x_2_dot, y_2_dot]


    def f(x_0, lmbda, T = 1):
        """Find f(x_0 + T)"""

        ### TODO: refactor, make into array, then transpose
        stepCnt = math.ceil(T / dt)

        # Need one more for the initial values
        ws = np.empty((stepCnt + 1, ))
        xs = np.empty((stepCnt + 1, ))
        ys = np.empty((stepCnt + 1, ))
        zs = np.empty((stepCnt + 1, ))

        # Setting initial values
        x_1 = x_0[0]
        y_1 = x_0[1]
        x_2 = x_0[2]
        y_2 = x_0[3]
        ws[0], xs[0], ys[0], zs[0] = x_1, y_1, x_2, y_2

        # Stepping through "time".
        for i in range(stepCnt):
            derivs = dots([ ws[i], xs[i], ys[i], zs[i] ], lmbda )

            ws[i + 1] = ws[i] + (derivs[0] * dt)
            xs[i + 1] = xs[i] + (derivs[1] * dt)
            ys[i + 1] = ys[i] + (derivs[2] * dt)
            zs[i + 1] = zs[i] + (derivs[3] * dt)

        return [ ws[-1], xs[-1], ys[-1], zs[-1] ]

    def g(x_0, lmbda, T = 1):
        """objective function"""
        return quad_sq_distance( f(x_0, lmbda, T), f(x_0, lmbda, 0) )

    def g_T(x_0):
        """g instantiated with a fixed period"""
        return g(x_0, lmbda, T)

    def newton_search(x_0, T = 1, N = 25):

        x = x_0

        hessian = nd.core.Hessian(g_T)
        jacobian = nd.core.Jacobian(g_T)

        for i in range(N):    
            adjust = np.matmul(np.linalg.inv(hessian(x)), np.transpose( jacobian(x)))
            adjust = np.transpose(adjust)[0]
            #print(x)
            #print(adjust)
            x = list_subtract(x, adjust)


        print(g_T(x))
        print(x)

    def plot_sim_path(x_0, T):
        stepCnt = math.ceil(T / dt)

        # Need one more for the initial values
        ws = np.empty((stepCnt + 1,))
        xs = np.empty((stepCnt + 1,))
        ys = np.empty((stepCnt + 1,))
        zs = np.empty((stepCnt + 1,))

        # Setting initial values
        x_1 = x_0[0]
        y_1 = x_0[1]
        x_2 = x_0[2]
        y_2 = x_0[3]
        ws[0], xs[0], ys[0], zs[0] = x_1, y_1, x_2, y_2

        # Stepping through "time".
        for i in range(stepCnt):
            # Derivatives of the W, X, Y, Z state
            derivs = dots([ ws[i], xs[i], ys[i], zs[i] ], lmbda )

            ws[i + 1] = ws[i] + (derivs[0] * dt)
            xs[i + 1] = xs[i] + (derivs[1] * dt)
            ys[i + 1] = ys[i] + (derivs[2] * dt)
            zs[i + 1] = zs[i] + (derivs[3] * dt)

        plot_quad(ws, xs, ys, zs, 0)
    
    if expmt == 'search':
        newton_search(start_pt)
    if expmt == 'plot':
        plot_sim_path(x_0, T)
    
experiment_1((10.2, 
              9.3, 
              14.4, 
              12.2) , expmt = 'plot')

experiment_1((4.2, 3.3, 4.4, 2.2),
             T = 10000, 
             lmbda = [0.086, 0.141, 0.773],
             expmt = 'plot')

experiment_1((4.2, 3.3, 4.4, 2.2),
             T = 1000, 
             lmbda = [0.086, 0.141, 0.773],
             expmt = 'search')



