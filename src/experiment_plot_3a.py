"""
experiment_plot_3a.py
plots from mintrick 3a
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

    evo = Evolution_ColluciNunez()

    print(evo)

    print(evo)
    print("============")

    expmt = Experiment( evo = evo, 
                        title = "Plot for Colucci Nunez solution 3a", 
                        descr = "For Colluci Nunez population dynamics system")

    start_pt = [6.5254581289417709, 0.24979777320847995, 0.75194036920525941, 4.2807722423952443]
    start_T = 40.0523045328004494

    # expmt.setParams(start_T = start_T, 
    #                 start_pt = start_pt )

    print("============")
    print(expmt)

    plt = evo.gen_plot( x_0 = start_pt, 
                        T = start_T,
                        stepCnt = 10000,
                        plot_type = 2, 
                        DEBUG = False)

    expmt.savePlot(plt)

    print(evo.F(x_0 = [6.4254581289417709, 0.23979777320847995, 0.74194036920525941, 4.2707722423952443], 
                T = 4.0523045328004494 ) )
    #6.7137523948875257, 0.25050456528021109, 0.7503460124735859, 4.4766520351302388, 3.9800410221931561
    plt.show()

if __name__=="__main__":
    main()
    
