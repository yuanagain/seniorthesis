## evolution.py
## Yuan Wang

from thesis_utils import *
from thesis_defaults import *
from thesis_poincare_utils import *
from thesis_plot_utils import *

class Evolution:
    """
    Evolution functions defining a dynamical systems on C^2. Need to specify a __call__ and __str__ method.
    """
    def __init__(self, dt = 0.01):
        self.dt = 0.01
        return

    def f(self, x_0):
        return self.__call__(x_0)

    def F(self, x_0, T = 1, dt = 0.01):
        """
        Returns path at time T starting from x_0, resolution dt.
        """

        if dt == None:
            dt = self.dt

        forwards = T >= 0 

        stepCnt = abs(math.ceil(T / self.dt))


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
            derivs = self.__call__([ ws[i], xs[i], ys[i], zs[i] ])

            if forwards:
                ws[i + 1] = ws[i] + (derivs[0] * dt)
                xs[i + 1] = xs[i] + (derivs[1] * dt)
                ys[i + 1] = ys[i] + (derivs[2] * dt)
                zs[i + 1] = zs[i] + (derivs[3] * dt)

            else: 
                ws[i + 1] = ws[i] - (derivs[0] * dt)
                xs[i + 1] = xs[i] - (derivs[1] * dt)
                ys[i + 1] = ys[i] - (derivs[2] * dt)
                zs[i + 1] = zs[i] - (derivs[3] * dt)

        return [ ws[-1], xs[-1], ys[-1], zs[-1] ]

    def plot_sim(self, x_0, T = 1000, plot_type = 0, dt = None):
        """
        Plots simulated path of system starting at x_0 for T, at resolution dt. 
        Plot types: see thesis_plot_utils.py
        """
        if dt == None:
            dt = self.dt

        plt = self.gen_plot(x_0 = x_0, T = T, plot_type = plot_type, dt = dt, )
        plt.show()

    def gen_plot(self, x_0, T = 1000.0, plot_type = 0, dt = None, stepCnt = None, DEBUG = False):
        """
        Generates a plot of simulated path of system starting at x_0 for T, at resolution dt. 
        Plot types: see thesis_plot_utils.py
        """
        if dt == None:
            dt = self.dt

        if stepCnt != None:
            dt = float(T) / stepCnt
            
        else:
            stepCnt = math.ceil(T / dt)

        if DEBUG:
            print("DEBUG: Step Count: " + str(stepCnt))

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
            derivs = self.__call__([ ws[i], xs[i], ys[i], zs[i] ] )

            ws[i + 1] = ws[i] + (derivs[0] * dt)
            xs[i + 1] = xs[i] + (derivs[1] * dt)
            ys[i + 1] = ys[i] + (derivs[2] * dt)
            zs[i + 1] = zs[i] + (derivs[3] * dt)

            if DEBUG:
                print( "DEBUG: Pt: "+ str([ws[i+1], xs[i+1], ys[i+1], zs[i+1]] ))

        txt = self.__str__() + '\n' + param_string(x_0, T, dt)
        
        return gen_plot(ws, xs, ys, zs, plot_type, txt = txt)



class Evolution_ColluciNunez(Evolution):

    def __init__(self,  e = 1.0, 
                        u_c = 1.0, 
                        u_g = 1.0, 
                        a_c = 3.0, 
                        a_g = 1.0, 
                        m = 1.0, 
                        dt = 0.01):

        self.dt = 0.01
        self.e, self.u_c, self.u_g, self.a_c, self.a_g, self.m = e, u_c, u_g, a_c, a_g, m
        return

    def __call__(self, x_0):
        e, u_c, u_g, a_c, a_g, m = self.e, self.u_c, self.u_g, self.a_c, self.a_g, self.m

        Z = x_0[0]
        C = x_0[1]
        G = x_0[2]
        P = x_0[3]

        Z_dot = Z * (-e + u_c * C + u_g * G)
        C_dot = C * (a_c * P - u_c * Z) - m * C * Z
        G_dot = G * (a_g * P - u_g * Z) + m * C * Z
        P_dot = P * (-a_g * G - a_c * C) + e * Z

        return [Z_dot, C_dot, G_dot, P_dot]

    # def f(self, x_0, T = 1, dt = 0.01):
    #     return super(Evolution_1a, self).f(x_0, T = 1, dt = 0.01)

    def __str__(self):

        return  "Z' = Z * (-e + u_c * C + u_g * G) \n" + \
                "C' = C * (a_c * P - u_c * Z) - m * C * Z \n" + \
                "G' = G * (a_g * P - u_g * Z) + m * C * Z \n" + \
                "P' = P * (-a_g * G - a_c * C) + e * Z \n" + \
                "e = " + str(self.e) + "; u_c = " +  str(self.u_c) + \
                "; u_g = " +  str(self.u_g) + "; a_c = " + \
                str(self.a_c) + "; a_g = " +  str(self.a_g) + \
                "; m = " + str(self.m)


class Evolution_Toy_A(Evolution):

    def __call__(self, x_0):

        x_1 = x_0[0]
        x_2 = x_0[1]
        x_3 = x_0[2]
        x_4 = x_0[3]

        x_1_dot = -x_1 + x_1 * x_2 + x_1 * x_3
        y_1_dot = 3 * x_2 * x_4 - 2 * x_1 * x_2
        x_2_dot = x_3 * x_4 - x_3 * x_1 + x_2 * x_1
        y_2_dot = -x_3 * x_4 - 3 * x_2 * x_4 + x_1

        return [x_1_dot, y_1_dot, x_2_dot, y_2_dot]

    # def f(self, x_0, T = 1, dt = 0.01):
    #     return super(Evolution_1a, self).f(x_0, T = 1, dt = 0.01)

    def __str__(self):

        return  "x1' = -x1 + x1*x2 + x1*x3 \n" + \
                "x2' = 3*x2*x4 - 2*x1*x2 \n" + \
                "x3' = x3*x4 - x3*x1 + x2*x1 \n" + \
                "x4' = -x3*x4 - 3*x2*x4 + x1 \n"



class Evolution_1a(Evolution):
    def __init__(self, lmbda, dt = 0.01):
        self.lmbda = lmbda
        self.dt = dt

    def __call__(self, x_0):

        x_1 = x_0[0]
        y_1 = x_0[1]
        x_2 = x_0[2]
        y_2 = x_0[3]

        x_1_dot = self.lmbda[1] * (x_1**2 - y_1**2) - (self.lmbda[1] + self.lmbda[2]) * (x_1*x_2 - y_1*y_2)
        y_1_dot = 2 * self.lmbda[1] * x_1 * y_1 - (self.lmbda[1] + self.lmbda[2]) * (x_1*y_2 + y_1*x_2)
        x_2_dot = self.lmbda[0] * (x_2**2 - y_2**2) - (self.lmbda[0] + self.lmbda[2]) * (x_1*x_2 - y_1*y_2)
        y_2_dot = 2 * self.lmbda[0] * x_2 * y_2 - (self.lmbda[0] +self.lmbda[2]) * (x_1*y_2 + y_1*x_2)

        return [x_1_dot, y_1_dot, x_2_dot, y_2_dot]

    # def f(self, x_0, T = 1, dt = 0.01):
    #     return super(Evolution_1a, self).f(x_0, T = 1, dt = 0.01)

    def __str__(self):
        return  "dz1/dt = lambda_2 * z1^2 - (lambda_2 + lambda_3) * z1 * z2 \n" + \
                "dz2/dt = lambda_1 * z2^2 - (lambda_1 + lambda_3) * z1 * z2 \n" + \
                "lambda_1: " + str(self.lmbda[0]) + \
                "; lambda_2: " + str(self.lmbda[1]) + "; lambda_3: " + str(self.lmbda[2])


class Evolution_Toy_B(Evolution):

    def __call__(self, x_0, mu = 3):
        """
        mu < 4
        """
        x_1 = x_0[0]
        x_2 = x_0[1]
        x_3 = x_0[2]
        x_4 = x_0[3]

        x_1_dot = x_2
        y_1_dot = mu * (1 - x_1**2) * x_2 - x_1
        x_2_dot = 0
        y_2_dot = 0

        return [x_1_dot, y_1_dot, x_2_dot, y_2_dot]

    # def f(self, x_0, T = 1, dt = 0.01):
    #     return super(Evolution_1a, self).f(x_0, T = 1, dt = 0.01)

    def __str__(self):

        return  "x1' = -x1 + x1*x2 + x1*x3 \n" + \
                "x2' = 3*x2*x4 - 2*x1*x2 \n" 



def main():
    ##
    ev = Evolution_1a(lmbda = lmbda_set_1)
    print(ev)
    print(ev(default_start))
    print(ev.F(default_start, T = 1000))
    ev.plot_sim(x_0 = default_start, plot_type = 1)

    ## Test toy
    ev = Evolution_Toy_A(dt = 0.5)
    print(ev)
    print(ev([2,2,2,2]))
    print(ev.F([2,2,2,2], T = 1))
    ev.plot_sim(x_0 = default_start, plot_type = 0)


if __name__=="__main__":
    main()

