## evolutionint.py
## Yuan Wang

from thesis_utils import *
from thesis_defaults import *
from thesis_poincare_utils import *
from thesis_plot_utils import *

from interval import interval, inf, imath
from intervals.IntervalN import IntervalN, _scalar

class EvolutionInt:
    """
    Evolution functions defining a dynamical systems on C^2. Need to specify a __call__ and __str__ method.
    """
    def __init__(self, hk =  [_scalar(0.01) for i in range(1000)] ):
        self.hk = hk
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

