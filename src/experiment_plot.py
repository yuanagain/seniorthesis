"""
experiment.py
Defines Experiment interface for our numerical experiments
Author: Yuan Wang
"""

from thesis_utils import *
from thesis_defaults import *
from thesis_poincare_utils import *
from thesis_plot_utils import *

from evolution import *
from experiment import *

import os
import datetime


def main():
    """
    Testing
    """

    print("============")
    evo = Evolution_1a(lmbda = lmbda_set_1)
    print(evo)
    print("============")
    print("============")
    expmt = Experiment(evo = evo, title = "Billiards Ball System Plot", descr = "Plot billiards ball system")
    print("============")
    print(expmt)

    plt = evo.gen_plot(default_start, plot_type = 1, stepCnt = 100000, DEBUG = False)
    expmt.savePlot(plt)

    plt3 = evo.gen_plot(default_start, plot_type = 0, stepCnt = 100000, DEBUG = False)
    expmt.savePlot(plt3)

    expmt.saveText("Text appended to info file successfully")

if __name__=="__main__":
    main()
    
