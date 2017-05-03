"""
experiment_intervalsim.py
Interval simulation, as in first part of Zgliczynski paper
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

from evolutioninterval import *
from experiment import *

class ExperimentIntervalSim(Experiment):

    # def __init__(self, evo, title = "TITLE", descr = ""):
    #     self.title = title
    #     self.descr = descr
    #     self.evo = evo
    #     self.setup()
    #     self.fig_ct = 0

    def setParams(self, step_ct = 1, start_pt = default_start):
        self.params['step_ct'] = step_ct
        self.params['start_pt'] = start_pt
        self.saveParams()

    def run(self):
        """
        Runs the Experiment
        """
        data = self.evo.generate(   self.params['start_pt'], 
                                    h = None,  
                                    p_e = 2,
                                    stepCt = self.params['step_ct'] )

        self.print(data)


def main():
    """
    Testing
    """

    print("============")
    evo = Evolution_Valdez(lmbda = lmbda_set_1)
    print(evo)
    print("============")
    print("============")

    expmt = ExperimentIntervalSim(evo = evo, 
                                title = "Zgliczynski Rigorous Interval Simulation", 
                                descr = "Rigorous Interval Simulation w/ Lohner type algo per Zgliczynski paper")

    expmt.setParams(step_ct = 100, start_pt = toArray(tupleToIntervalVector(default_start)) )

    # print("============")
    expmt.run()

    

if __name__=="__main__":
    main()