"""
experiment_mintrick2.py
Defines Experiment interface for our numerical experiments
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

from evolution import *
from experiment import *

class ExperimentMinTrick2(Experiment):

    # def __init__(self, evo, title = "TITLE", descr = ""):
    #     self.title = title
    #     self.descr = descr
    #     self.evo = evo
    #     self.setup()
    #     self.fig_ct = 0

    def setParams(self, T = 1, start_pt = default_start):
        self.params['T'] = T
        self.params['start_pt'] = start_pt
        self.saveParams()

    def run(self):
        """
        Runs the Experiment
        """
        self.print("Running Newton Search, varying t")
        self.newton_search()

    def f_integrand(self, t = 0):
        """
        Distance |f_T+t(x_0) - f_t(x_0)|^2, = 0 iff orbit of period T exists starting at f_t(x_0)
        """
        return quad_sq_distance(self.evo.F(self.params['start_pt'], t), 
                                self.evo.F(self.params['start_pt'], t + self.params['T']))
    
    def phi(self, t):
        """
        Integrated |f_T+t(x_0) - f_t(x_0)|^2, = 0 iff orbit of period T exists starting at f_t(x_0)
        """
        return quad(lambda t: self.f_integrand(t), 0, self.params['T'])[0]
    
    def newton_search(self, t = None):
        """
        Minimize f_integrand while varying t
        """
        if t == None:
            t = self.params['T'] / 2


        self.print("min_t:")
        min_t = newton(self.phi, t)
        self.print(min_t)
        self.print("F(min_t)")
        self.print(self.evo.F(self.params['start_pt'], min_t))

def main():
    """
    Testing
    """

    print("============")
    #evo = Evolution_1a(lmbda = lmbda_set_1)
    evo = Evolution_Toy_A()
    print(evo)
    print("============")
    print("============")

    expmt = ExperimentMinTrick2(evo = evo, 
                                title = "Minimization Trick 1", 
                                descr = "Seeking orbits by minimizing |f_T(x) - x|^2")

    expmt.setParams(T = 1, start_pt = default_start)

    print("============")
    print(expmt)

    plt = evo.gen_plot(default_start, plot_type = 0)
    expmt.savePlot(plt)
    expmt.run()

    

if __name__=="__main__":
    main()