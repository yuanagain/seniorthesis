"""
experiment_c1lohner.py
Poincare C^1 Lohner Algorithm
Author: Yuan Wang

On Lohner Algorithm
https://books.google.com/books?id=7a-8yyjQVLcC&pg=PA178&lpg=PA178&dq

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

class ExperimentC1Lohner(Experiment):

    def setParams(self, T = 1000, start_pt = default_start):
        self.params['T'] = T
        self.params['start_pt'] = start_pt
        self.saveParams()

    def run(self, T = None, dt = 0.01, stepCnt = 10000):



def main():
    """
    Testing
    """

    print("============")
    #evo = Evolution_1a(lmbda = lmbda_set_1)
    evo = Evolution_ColluciNunez()

    print(evo)

    expmt = ExperimentC1Lohner( evo = evo, 
                                title = "Poincare C^1 Lohner Algorithm", 
                                descr = "C^1 Lohner algorithm for period detection")

    # expmt.setParams(T = 4, start_pt = default_start)

    expmt.setParams(hyperplane = HyperPlane(-4, 12, 2, -10, -1.2), 
                    T = 30, 
                    start_pt = [9.1, 4.1, 3.2, 4.5] )

    print("============")
    print(expmt)

    expmt.run(T = None, stepCnt = 1000)

if __name__=="__main__":
    main()

