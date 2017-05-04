"""
experiment_poincare.py
Poincare map generation
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

class ExperimentPoincare(Experiment):

    def setParams(self, hyperplane = HyperPlane(1, 1, 1, 1, 4), T = 1000, start_pt = default_start):
        self.hyperplane = hyperplane
        self.params['T'] = T
        self.params['start_pt'] = start_pt
        
        self.saveParams()

    def run(self, T = None, dt = 0.01, stepCnt = 10000):
        """Simulate path, collecting Poincare crossings"""

        start_pt = self.params['start_pt']

        if T == None:
            T = self.params['T']

        if stepCnt != None:
            dt = float(T) / stepCnt
            
        else:
            stepCnt = math.ceil(T / dt)


        # dt = 0.01
        # stepCnt = 100000

        # Need one more for the initial values
        ws = np.empty((stepCnt + 1,))
        xs = np.empty((stepCnt + 1,))
        ys = np.empty((stepCnt + 1,))
        zs = np.empty((stepCnt + 1,))

        crossings = np.empty((stepCnt + 1,))
        pts = np.empty((stepCnt + 1,))

        ws[0], xs[0], ys[0], zs[0] = start_pt[0], start_pt[1], start_pt[2], start_pt[3]

        current_pt = list(start_pt)
        
        crossings[0] = 0

        pts[0] = self.hyperplane(current_pt)

        intersect_checker = IntersectChecker(self.hyperplane)
        trace = [ws, xs, ys, zs]

        ## for tracking min/max/mean of path, relative to hyperplane
        
        # Stepping through "time".
        self.print("\n\nCrossings:")
        for i in range(stepCnt):
            # Derivatives of the W, X, Y, Z state

            derivs = self.evo(current_pt)

            old_pt = list(current_pt)
            ## compute new point
            for j in range(4):
                trace[j][i + 1] = old_pt[j] + (derivs[j] * dt)
                current_pt[j] = trace[j][i + 1]

            pts[i + 1] = self.hyperplane(current_pt)

            crossings[i + 1] = intersect_checker(current_pt)
            # print(hyperplane(pt))

            if crossings[i + 1] != 0:
                self.print((ws[i + 1], xs[i + 1], ys[i + 1], zs[i + 1]))

        self.print("\nMax:")
        self.print(max(pts))
        self.print("Min:")
        self.print(min(pts))
        self.print("Av:")
        self.print(sum(pts) / len(pts))

        ws, xs, ys, zs = poincareExtract(ws, xs, ys, zs, crossings)
        
        # for i in range(len(ws)):
        #     self.print( "(" + str(ws[i]) + ", " + str(xs[i]) + ", " + str(ys[i]) + ", " + str(zs[i]) + ")" )

        self.savePlot(poincarePlot(ws, xs, ys, zs, str(self.hyperplane)))

        # if expmt == 'accumulate':
        #     return [ws, xs, ys, zs, crossings]
        

def main():
    """
    Testing
    """

    print("============")
    #evo = Evolution_1a(lmbda = lmbda_set_1)
    evo = Evolution_ColluciNunez()

    print(evo)

    expmt = ExperimentPoincare( evo = evo, 
                                title = "Poincare map generation", 
                                descr = "Leveraging Poincare maps to gain insights about our system")

    # expmt.setParams(T = 4, start_pt = default_start)

    expmt.setParams(hyperplane = HyperPlane(-4, 12, 2, -10, -1.2), 
                    T = 30, 
                    start_pt = [9.1, 4.1, 3.2, 4.5] )

    print("============")
    print(expmt)

    expmt.run(T = None, stepCnt = 1000)

if __name__=="__main__":
    main()

