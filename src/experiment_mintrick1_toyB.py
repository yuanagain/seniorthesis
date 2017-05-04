"""
experiment_mintrick1_toyB.py
Minimization trick 1, toy problem B
Author: Yuan Wang
"""

from thesis_utils import *
from thesis_defaults import *
from thesis_poincare_utils import *
from thesis_plot_utils import *

import scipy.integrate as integrate
import scipy.special as special
from scipy.integrate import quad
from scipy.optimize import newton

import numdifftools as nd

from evolution import *
from experiment import *

class ExperimentMinTrick1(Experiment):

    def setParams(self, T = 1, start_pt = default_start):
        self.params['T'] = T
        self.params['start_pt'] = start_pt
        self.saveParams()

    def run(self):
        """
        Runs the Experiment
        """
        self.print("Running Newton Search, Varying x_0")
        self.newton_search()

    def g(self, x_0):
        """
        Function we are trying to minimize
        """
        return quad_sq_distance(self.evo.F(x_0, self.params['T']), 
                                self.evo.F(x_0, 0) )

    def phi(self, t):
        """
        Integrated |f_T+t(x_0) - f_t(x_0)|^2, = 0 iff orbit of period T exists starting at f_t(x_0)
        """
        return quad(lambda t: self.f_integrand(t), 0, self.params['T'])[0]
    
    def newton_search(self, t = None):
        """
        Minimize f_integrand while varying x_0
        """
        hessian = nd.core.Hessian(self.g)
        jacobian = nd.core.Jacobian(self.g)

        x = self.params['start_pt']

        self.print("*(x, g(x)):*")
        for i in range(50):    
            adjust = np.matmul(np.linalg.inv(hessian(x)), np.transpose( jacobian(x)))
            adjust = np.transpose(adjust)[0]

            #print(x)
            #print(adjust)

            # self.print("SUM(X)=" + str(sum(x)))
            x = list_subtract(x, adjust)

            self.print((x, self.g(x)))
        

def main():
    """
    Testing
    """

    print("============")
    #evo = Evolution_1a(lmbda = lmbda_set_1)
    evo = Evolution_Toy_B()

    print(evo)
    print("============")
    print("============")

    expmt = ExperimentMinTrick1(evo = evo, 
                                title = "Minimization Trick 1B", 
                                descr = "Seeking orbits by minimizing |f_T(x) - x|^2")

    expmt.setParams(T = 2, start_pt = [1,2,0,0])
    print("============")
    print(expmt)

    plt = evo.gen_plot( x_0 = [1,2,0,0], 
                        T = 1,
                        dt = 0.1,
                        plot_type = 0)
    expmt.savePlot(plt)


    expmt.run()



if __name__=="__main__":
    main()