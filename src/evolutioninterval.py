## evolutioninterval.py
## Yuan Wang

from thesis_utils import *
from thesis_defaults import *
from thesis_poincare_utils import *
from thesis_plot_utils import *
from thesis_matrix_utils import *

import time

import numpy as np
import math

from intervals.intervaln import *
grow_int = interval([1.0 - 1e-7, 1.0 + 1e-7])

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

        x_1 = interval(x_0[0][0])
        y_1 = interval(x_0[1][0])
        x_2 = interval(x_0[2][0])
        y_2 = interval(x_0[3][0])

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


    def refine(self, x_k, h_k, Y_0):
        """
        Refines enclosure of [x_k] to W_1 satisfying condition (22)
        """
        Y_i = Y_0
        h_k_int = interval([0, h_k])

        def condition(x_k, h_k, Y_i):
            return contains(Y_i, x_k + self.__call__(Y_i) * h_k_int )

        while not condition(x_k, h_k, Y_i):

            Y_i = x_k + h_k_int * self.__call__(Y_i)

        return Y_i

    def nabla(self, x):
        """
        Returns jac(x)
        """
        # x = [y[0] for y in x]
        x_1, x_2, x_3, x_4 = None, None, None, None

        ### Matrix indexing annoyingly incompatible with double array
        if type(x) == np.matrix:
            x_1 = interval(x[0, 0])
            x_2 = interval(x[1, 0])
            x_3 = interval(x[2, 0])
            x_4 = interval(x[3, 0])

        else:
            x_1 = interval(x[0][0])
            x_2 = interval(x[1][0])
            x_3 = interval(x[2][0])
            x_4 = interval(x[3][0])

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

        return toArray(jac)


    def f(self, x_0):
        ## column matrix
        return self.__call__(x_0)
        #return toArray([ [intvl] for intvl in self.__call__(x_0) ])

    def F(self, x_0, h = None, N = 10):
        """
        Returns path at time T starting from x_0, resolution dt.
        """

        ## uniform spacing by default
        if h == None:
            h = [self.dt for i in range(N) ]

        # Need one more for the initial values

        # Setting initial values
        # print("======")
        xpts = np.empty(shape=(N + 1,), dtype = interval)
        # print("Xpts:")
        # print(xpts)
        xpts[0] = (toArray(x_0))
        # Stepping through "time".
        for i in range(N):
            h_k = h[i]
            x_k = xpts[i]

            # print("x_k = "  + str(x_k))
            # print("TYPE x_k: " + str(type(x_k)) )

            derivs = self.__call__(x_k)
            scaled = derivs * h_k
            #next = safeAdd(x_k, scaled)
            # print("SHAPE scaled: " + str(scaled.shape))
            # print("SHAPE x_k: " + str(x_k.shape))
            # #next = np.matrix(x_k) + np.matrix(scaled)
            next = x_k + scaled
            # print('NEXT: ' + str(type(next)))

            xpts[i+1] = next
            # print("TYPE x_inserted: " + str(type(xpts[i+1])) )

            ### Get W_1 via refinement

            ### Evaluate taylor approx

            
            ### Compute solution for next time step

        return xpts

    def generate(self, x_0, h = None, p_e = 2, stepCt = 1000):
        """
        generate rigorous bounds of solution [x_k] via C_0 Lohner algo
        """

        ## uniform spacing by default
        if h == None:
            h = [self.dt for i in range(stepCt) ]

        # Need one more for the initial values

        # Setting initial values
        xpts = np.empty(shape=(stepCt + 1,), dtype = interval)
        xpts[0] = np.matrix((toArray(x_0)))
        # Stepping through "time".
        for i in range(stepCt):
            h_k = h[i]
            x_k = toArray(xpts[i])

            ### x_k of the form [[interval],[interval],[interval],[interval]]

            # derivs = self.__call__(x_k)
            # scaled = derivs * h_k
            # #next = safeAdd(x_k, scaled)
            # next = np.matrix(x_k) + np.matrix(scaled)
            # xpts = xpts + [next]

            ### Get W_1 via refinement
            grow = interval([1.0 - 1e-7, 1.0 + 1e-7])
            Y_0 = [[item[0] * grow_int] for item in x_k]
            W_1 = self.refine(x_k, h_k, Y_0)

            ### Evaluate taylor approx

            ## compute A_k, an interval matrix that is the Jacobian of Phi evaluated on [x_k].
            A_k = self.nablaPhi(h = h_k, x = x_k, p = p_e)

            ### Compute solution for next time step
            x_k_next =  self.Phi(h = h_k, x = midpoint(x_k), p = p_e) + \
                        self.dfdt(W_1, p = p_e) * (h_k**(p_e + 1) / math.factorial(p_e + 1)) + \
                        A_k * (x_k - midpoint(x_k))

            ### Store and repeat
            xpts[i+1] = x_k_next
            # break
            

        return xpts

    def dfdt(self, x, p = 2):
        """
        derivative of evolution function with respect to time at x of order 2
        """
        return self.dphidt(t = 0, x = x, p = p + 1)

    def phi(self, t, x_0):
        """
        Just an alias for evolution
        """
        if t == 0:
            return x_0

        else:
            raise ValueError("We haven't solved phi yet!")

    def dphidt(self, t, x, p = 1):
        """
        Returns the pth order derivative of phi with respect to time. Returns as matrix
        """
        if p == 0:
            return np.matrix( self.phi(t = t, x_0 = x) )

        if p == 1:
            if t == 0:
                return np.matrix( self.f(x) )
            else:
                raise ValueError("We haven't solved phi yet!")

        if p == 2:
            if t == 0:
                ## though self.phi(t = t, x_0 = x) = x, we're keeping it in this form to 
                ## make it more interpretable
                return np.matrix( self.nabla( self.phi(t = t, x_0 = x) ) ) *  \
                        np.matrix( self.f( self.phi(t = t, x_0 = x) ) )

            else:
                raise ValueError("We haven't solved phi yet!")

        if p == 3:
            if t == 0:
                return  np.matrix(self.nabla(self.dphidt(t, x, p = 1))) * \
                            np.matrix(self.f(self.phi(t, x))) + \
                        np.matrix(self.nabla(self.phi(t, x))) * \
                            np.matrix(self.nabla(self.phi(t, x))) * \
                            np.matrix(self.f(self.phi(t, x)))

            else:
                raise ValueError("We haven't solved phi yet!")

        else:
            raise ValueError("Deriviative of phi of order p = " + str(p) + " unimplemented")

    def Phi(self, h, x, p = 2):
        """
        Taylor approximation of order p at x, time step h
        """

        summand = _zero()
        for i in range(p + 1):
            summand = summand + (h**i / math.factorial(i)) * self.dphidt(t = 0, x = x, p = i)

        return summand

    def V(t, x):
        jac = self.nabla(x)
        return np.matrix(jac) * np.matrix(self.phi(t, x))

    def dVdt(self, t, x, p):
        """
        Jacobian matrix of phi with respect to initial state x_0.
        """
        if p == 0:
            return np.matrix(_oneArr())

        elif p == 1:
            if t == 0:
                jac = self.nabla(x)
                return np.matrix(jac)

            else:
                raise ValueError("We haven't solved phi yet!")

        elif p == 2: 
            if t == 0:
                jac = np.matrix(self.nabla(x))
                return np.matrix( self.nabla(np.matrix(self.f(x))) ) *  + \
                        jac * jac

            else:
                raise ValueError("We haven't solved phi yet!")

        else:
            raise ValueError("Deriviative of V of order p = " + str(p) + " unimplemented")

    def nablaPhi(self, h, x, p = 2):
        """
        Jacobian of Phi evaluated at x
        """
        summand = _zeroArr()
        for i in range(p + 1):
            summand = summand + (h**i / math.factorial(i)) * self.dVdt(t = 0, x = x, p = i)

        return summand

def _oneArr():
    return toArray([[1,0,0,0],
                    [0,1,0,0],
                    [0,0,1,0],
                    [0,0,0,1]])

def _zeroArr():
    return toArray([[0,0,0,0],
                    [0,0,0,0],
                    [0,0,0,0],
                    [0,0,0,0]])

def _zero():
    return toArray(tupleToIntervalVector([0,0,0,0]))

def main():
    ##
    ev = Evolution_Valdez(lmbda = lmbda_set_1)
    print(ev)
    print(ev([interval(pt) for pt in default_start]))

    intvl_start = listToCol([interval(pt) for pt in default_start])
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

    

    Y_0 = [[item[0] * grow_int] for item in x_k]
    print("Y_0 = ")
    print(Y_0)
    # Y_0 = grow * x_k
    # print("Y_0 = ")
    # print(Y_0)
    print("REFINING")
    print(ev.refine(x_k, h_k, Y_0))

    print(( np.matrix(jac_x) * np.matrix(start) ).shape)

    print("testing midpoint")
    print(midpoint(Y_0))

    print("Testing dphi/dt")
    for p_e in range(4):
        print("p_e = " + str(p_e))
        testdphidt = ev.dphidt(t = 0, x = x_k, p = p_e)
        print("VALIDATED?: " + str(validateType(testdphidt, show = True)))
        # print(testdphidt)
        # print(type(testdphidt))
        # print(type(testdphidt[0]))
        # print(type(testdphidt[0, 0]))

    print("Testing df/dt")
    for p_e in range(3):
        print("p_e = " + str(p_e))
        testdfdt = ev.dfdt(x_k, p = p_e)
        print("VALIDATED?: " + str(validateType(testdfdt, show = True)))
        # print(testdfdt)
        # print(type(testdfdt))
        # print(type(testdfdt[0]))
        # print(type(testdfdt[0, 0]))

    print("Testing dV/dt")
    for p in range(3):
        print("p_e = " + str(p))
        testdvdt = ev.dVdt(t = 0, x = x_k, p = p)
        print("VALIDATED?: " + str(validateType(testdvdt, show = True)))
        # print(type(testdvdt))
        # print(type(testdvdt[0]))
        # print(type(testdvdt[0, 0]))

    print("Testing Phi")
    testPhi = ev.Phi(h = 0.01, x = x_k, p = 2)
    validateType(testPhi, show = True)

    print("Testing nablaPhi")
    testnablaPhi = ev.nablaPhi(h = 0.01, x = x_k, p = 2)
    print("VALIDATED?: " + str(validateType(testnablaPhi, show = True)))


if __name__=="__main__":
    main()

