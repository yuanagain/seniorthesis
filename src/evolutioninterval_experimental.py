## evolutioninterval.py
## Yuan Wang

from thesis_utils import *
from thesis_defaults import *
from thesis_poincare_utils import *
from thesis_plot_utils import *
from thesis_matrix_utils import *

import numpy as np

from intervals.intervaln import *

class EvolutionInterval:
    """
    Evolution functions defining a dynamical systems on C^2. Need to specify a __call__ and __str__ method.
    """
    def __init__(self, dt = 0.01):
        self.dt = 0.01
        return

    def f(self, x_0):
        return self.__call__(x_0)

    def F(self, x_0, h = None, dt = 0.01):
        """
        Returns path at time T starting from x_0, resolution dt.
        """

        ## uniform spacing by default
        if h == None:
            h = [self.dt for i in range(10) ]

        # Need one more for the initial values

        # Setting initial values
        print("======")
        xpts = [x_0]
        # Stepping through "time".
        for h_k in h:
            derivs = self.__call__(xpts[-1])
            xpts = xpts + [xpts[-1] + derivs * h_k]

        return xpts


class Evolution_Valdez(EvolutionInterval):
    def __init__(self, lmbda, dt = 0.01):
        self.lmbda = lmbda
        self.dt = dt

    def __call__(self, x_0):
        # print("x_0 x_0 x_0 x_0 x_0 x_0 x_0 x_0 x_0 x_0 ")
        # print(x_0)

        x_1 = interval(x_0.x[0][0])
        y_1 = interval(x_0.x[1][0])
        x_2 = interval(x_0.x[2][0])
        y_2 = interval(x_0.x[3][0])

        x_1_dot = self.lmbda[1] * (x_1**2 - y_1**2) - (self.lmbda[1] + self.lmbda[2]) * (x_1*x_2 - y_1*y_2)
        y_1_dot = 2 * self.lmbda[1] * x_1 * y_1 - (self.lmbda[1] + self.lmbda[2]) * (x_1*y_2 + y_1*x_2)
        x_2_dot = self.lmbda[0] * (x_2**2 - y_2**2) - (self.lmbda[0] + self.lmbda[2]) * (x_1*x_2 - y_1*y_2)
        y_2_dot = 2 * self.lmbda[0] * x_2 * y_2 - (self.lmbda[0] +self.lmbda[2]) * (x_1*y_2 + y_1*x_2)

        res = toArray(listToCol( [x_1_dot, y_1_dot, x_2_dot, y_2_dot] ))

        return res

    # def f(self, x_0, T = 1, dt = 0.01):
    #     return super(Evolution_1a, self).f(x_0, T = 1, dt = 0.01)

    def __str__(self):
        return  "dz1/dt = lambda_2 * z1^2 - (lambda_2 + lambda_3) * z1 * z2 \n" + \
                "dz2/dt = lambda_1 * z2^2 - (lambda_1 + lambda_3) * z1 * z2 \n" + \
                "lambda_1: " + str(self.lmbda[0]) + \
                "; lambda_2: " + str(self.lmbda[1]) + "; lambda_3: " + str(self.lmbda[2])

    def nabla(self, x):
        """
        Returns jac(x)
        """
        # x = [y[0] for y in x]
        x_1 = interval(x.x[0][0])
        x_2 = interval(x.x[1][0])
        x_3 = interval(x.x[2][0])
        x_4 = interval(x.x[3][0])

        p_1 = self.lmbda[1]
        p_2 = self.lmbda[1] + self.lmbda[2]
        p_3 = self.lmbda[0]
        p_4 = self.lmbda[0] + self.lmbda[2]

        jac = [ [2*p_1*x_1 - p_2*x_3, -2*p_1*x_2 - p_2*x_4, -p_2*x_1, p_2*x_2 ] ,
                [2*p_2*x_1 - p_2*x_4, -2*p_1*x_1 - p_2*x_3, -p_2*x_2, -p_2*x_1 ] ,
                [-p_4*x_3, -p_4*x_4, 2*p_3*x_3 - p_4*x_1, -2*p_3*x_4 - p_4*x_2 ] ,
                [-p_4*x_4,  p_4*x_3, 2*p_3*x_4 - p_4*x_2,  2*p_3*x_3 - p_4*x_1 ] 
              ]

        arr = np.empty(shape=(4,4), dtype=interval)

        for i in range(len(arr)):
            for j in range(len(arr)):
                arr[i, j] = jac[i][j]

        res = toArray(jac)

        return res


    def f(self, x_0):
        ## column matrix
        return toArray([ [intvl] for intvl in self.__call__(x_0).x ])

    def F(self, x_0, h = None):
        """
        Returns path at time T starting from x_0, resolution dt.
        """

        ## uniform spacing by default
        if h == None:
            h = [self.dt for i in range(10) ]

        # Need one more for the initial values

        # Setting initial values
        print("======")
        xpts = [x_0]
        # Stepping through "time".
        for h_k in h:
            x_k = xpts[-1]

            derivs = self.__call__(x_k)
            scaled = derivs * h_k
            next = safeAdd(x_k, scaled)
            xpts = xpts + [next]

            ### Get W_1 via refinement

            ### Evaluate taylor approx

            
            ### Compute solution for next time step

        return xpts

    def generate(self, x_0, h = None):
        """
        generate solution via C_0 Lohner algo
        """

        ## uniform spacing by default
        if h == None:
            h = [self.dt for i in range(10) ]

        # Need one more for the initial values

        # Setting initial values
        print("======")
        xpts = [x_0]
        # Stepping through "time".
        for h_k in h:
            x_k = xpts[-1]

            derivs = self.__call__(x_k)
            scaled = derivs * h_k

            next = safeAdd(x_k, scaled)
            xpts = xpts + [next]

            ### Get W_1 via refinement



            ### Evaluate taylor approx

            
            ### Compute solution for next time step

        return xpts

    def Phi(self, h, z_k, p):

        return 

    def refine(self, x_k, h_k, Y_0):
        """
        Refines enclosure of [x_k] to W_1 satisfying condition (22)
        """
        Y_i = Y_0
        h_k_int = interval([0, h_k])

        def condition(x_k, h_k, Y_i):
            intvln_1 = IntervalN(colToRow(Y_i))
            intvln_2 = IntervalN(colToRow( x_k + h_k_int * self.f(Y_i) ))
            return intvln_1 in intvln_2

        while not condition(x_k, h_k, Y_i):

            Y_i = x_k + h_k_int * self.f(Y_i)

        return Y_current

def main():
    ##
    ev = Evolution_Valdez(lmbda = lmbda_set_1)
    print(ev)
    print(ev([interval(pt) for pt in default_start]))

    intvl_start = IntervalN(default_start)
    print("interval start")
    print(intvl_start)
    jac_x = ev.nabla(intvl_start)
    print('interval jacobian')
    print(jac_x)
    print('interval jacobian [0, ]')
    print(jac_x[0])
    print('interval jacobian [0,0]')
    print(jac_x[0][0])
    print(type(jac_x[0][0]))

    print("------------------")
    start = toArray(tupleToIntervalVector(default_start))
    print(start)
    print("ev(start)")
    print(ev(start))
    deriv = ev(start)
    scaled = deriv * 0.01
    print(type(start))
    print(type(scaled))
    sum = scaled + start
    print(sum)
    print(type(sum))

    print("F F F F F F F F F F F F F F F ")

    h_k = 0.01
    x_k = start

    print(type(x_k[0][0]))

    derivs = ev(x_k)
    print('derivs derivs derivs derivs')
    print(type(derivs))
    print(derivs)

    scaled = derivs * h_k
    print("scaled scaled scaled scaled")
    print(type(scaled))
    print(scaled)

    next = x_k + scaled
    print("next next next next")
    print(type(next))
    print(next)
    

    lst = ev.F(tupleToIntervalVector(default_start) )
    print('lst lst lst lst lst')
    print(len(lst))
    print(lst[2])
    print(type(lst))
    print(type(lst[2]))
    print(type(lst[2][3]))
    print(type(lst[2][3][0]))

    grow = interval([0.999, 1.001])

    Y_0 = [item[0] * grow for item in x_k]
    print("Y_0 = ")
    print(Y_0)
    Y_0 = grow * x_k
    print("Y_0 = ")
    print(Y_0)
    ev.refine(x_k, h_k, Y_0)


if __name__=="__main__":
    main()

