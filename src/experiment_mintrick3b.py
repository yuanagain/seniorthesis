"""
experiment_mintrick3b.py
Minimization trick 3: applied to Colucci, Nunez paper
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

class ExperimentMinTrick3(Experiment):

    def setParams(self, start_T = 1, start_pt = default_start):
        self.params['start_T'] = start_T
        self.params['start_pt'] = start_pt
        self.saveParams()

    def run(self):
        """
        Runs the Experiment
        """
        self.print("Running Newton Search, Varying x_0")
        self.newton_search()

    def g(self, x):
        """
        Function we are trying to minimize, 
        As we are varying period length, we have a function of (x_1, x_2, x_3, x_4, T)
        """

        x_0 = x[0:4]
        T = x[4]

        return quad_sq_distance(self.evo.F(x_0, T), 
                                self.evo.F(x_0, 0) )
    
    def newton_search(self, t = None):
        """
        Minimize f_integrand while varying x_0
        """
        hessian = nd.core.Hessian(self.g)
        jacobian = nd.core.Jacobian(self.g)

        x = self.params['start_pt'] + [self.params['start_T']]

        self.print("*(x, g(x)):*")
        for i in range(10):    
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
    evo = Evolution_ColluciNunez()

    print(evo)

    print(evo)
    print("============")

    expmt = ExperimentMinTrick3(evo = evo, 
                                title = "Minimization Trick 3b", 
                                descr = "For Colluci Nunez population dynamics system")

    expmt.setParams(start_T = 1.02, start_pt = [13.812112161101782, 0.084080320330548491, 0.60914816966688357, 6.3946593489005625] )
    print("============")
    print(expmt)

    plt = evo.gen_plot( x_0 = [13.812112161101782, 0.084080320330548491, 0.60914816966688357, 6.3946593489005625], 
                        T = 1.02,
                        stepCnt = 100000,
                        plot_type = 1, 
                        DEBUG = False)

    expmt.savePlot(plt)
    plt.show()
    #expmt.run()
    print(evo.F(x_0 = [13.812112161101782, 0.084080320330548491, 0.60914816966688357, 6.3946593489005625], 
            T = 1 ) )
    print(expmt.g([13.812112161101782, 0.084080320330548491, 0.60914816966688357, 6.3946593489005625, 1.02]))
    #6.7137523948875257, 0.25050456528021109, 0.7503460124735859, 4.4766520351302388, 3.9800410221931561
    


if __name__=="__main__":
    main()