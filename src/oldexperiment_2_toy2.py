## experiment_2.py

import scipy.integrate as integrate
import scipy.special as special
from scipy.integrate import quad
from scipy.optimize import newton

experiment_2_DEBUG = False

def experiment_2(start_pt = default_start, 
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
    
        #return [-x_1_dot, -y_1_dot, -x_2_dot, -y_2_dot]


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

    def f_integrand(t, x_0, lmbda, T = 1):
        return quad_sq_distance(f(x_0, lmbda, t + T), f(x_0, lmbda, t))
    
    def phi(t, x_0, lmbda):
        """What we want to minimize"""
        return quad(f_integrand, 0, T, args=(x_0, lmbda))[0]
    
    def phi_instance(t):
        return phi(t, start_pt, lmbda)
    
    def newton_search(t, T = 1, N = 25):
        newton(phi_instance, t)

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
        newton_search(t = 100)
    if expmt == 'plot':
        plot_sim_path(x_0, T)


def main():
    

if __name__=="__main__":
    main()





