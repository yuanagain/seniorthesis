"""
experiment_timing.py
Timing
Author: Yuan Wang
"""

from thesis_utils import *
from thesis_defaults import *
from thesis_poincare_utils import *
from thesis_plot_utils import *

from evolution import *
from experiment import *


import time
import os
import datetime


def main():
    """
    Testing
    """

    print("============")
    evo = Evolution_1a(lmbda = lmbda_set_1)
   	
    print("Timing")
    start = time.time()
    evo.gen_plot(x_0 = default_start, T = 1.0, stepCnt = 1000000)
    end = time.time()
    print("Time:")
    print(end - start)

if __name__=="__main__":
    main()
   

    