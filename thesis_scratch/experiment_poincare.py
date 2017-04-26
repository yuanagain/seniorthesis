## Poincare experiment
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math

from thesis_utils import *
from thesis_defaults import *
from thesis_poincare_utils import *
from thesis_plot_utils import *

import numdifftools as nd

dt = 0.01
DEBUG = False


import random

def experiment3(start_pt = default_start, 
                 T = 10000, 
                 lmbda = [default_lambda_1, default_lambda_2, default_lambda_3], 
                 res = 0.001,
                 hyperplane = HyperPlane(0.2, 0.4, 1.2, -2.1, -1.1),
                 expmt = "accumulate"):
    
    ## define evaluation function

    
        ## if reversing time
        #return [-x_1_dot, -y_1_dot, -x_2_dot, -y_2_dot]

    def run(x_0, T):
        """Simulate path, collecting Poincare crossings"""
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

        print(stepCnt)
        print(dt)

        # # Setting initial values
        # x_1 = x_0[0]
        # y_1 = x_0[1]
        # x_2 = x_0[2]
        # y_2 = x_0[3]
        ws[0], xs[0], ys[0], zs[0] = start_pt[0], start_pt[1], start_pt[2], start_pt[3]
        current_pt = list(start_pt)
        crossings[0] = 0
        pts[0] = hyperplane(current_pt)

        intersect_checker = IntersectChecker(hyperplane)
        trace = [ws, xs, ys, zs]

        ## for tracking min/max/mean of path, relative to hyperplane
        
        # Stepping through "time".
        for i in range(stepCnt):
            # Derivatives of the W, X, Y, Z state

            derivs = dots(current_pt, lmbda )
            old_pt = list(current_pt)
            ## compute new point
            for j in range(4):
                trace[j][i + 1] = old_pt[j] + (derivs[j] * dt)
                current_pt[j] = trace[j][i + 1]

            pts[i + 1] = hyperplane(current_pt)

            
            crossings[i + 1] = intersect_checker(current_pt)
            # print(hyperplane(pt))

            if crossings[i + 1] != 0:
                print(i)
                print((ws[i + 1], xs[i + 1], ys[i + 1], zs[i + 1]))

        print(max(pts))
        print(min(pts))
        print(sum(pts) / len(pts))

        poincareExtract(ws, xs, ys, zs, crossings)

        if expmt == 'trace':
            plot_quad(ws, xs, ys, zs, 1, txt = "")
  
        
        if expmt == 'plot':
            poincarePlot(ws, xs, ys, zs, crossings, str(hyperplane))

        if expmt == 'accumulate':
            return [ws, xs, ys, zs, crossings]

 
    return run(start_pt, T)

# experiment_3((0.032, 0.308, -0.1, -0.5),
#              T = 1083, 
#              lmbda = [0.086, 0.141, 0.773],
#              hyperplane = HyperPlane(4, -3, -1, -4, 0),
#              expmt = 'plot')

import random

# experiment_3(random_start,
#              T = 1083, 
#              lmbda = [0.086, 0.141, 0.773],
#              hyperplane = HyperPlane(1, 1, 1, 1, 1),
#              expmt = 'accumulate')

def experiment4():
    """
    Accumulate hyperplane stuff
    """
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    plt.subplots_adjust(top=0.85)
    
    # plt.xlim( (-10, 10) )
    # plt.ylim( (-10, 10) )
    # # plt.zlim( (-10, 10) )
    ax.set_xlim(-3, 3 )
    ax.set_ylim(-3, 3)
    ax.set_zlim(-3, 3)
    ## execute

    hyperplane = HyperPlane(1, 1, 1, 1, 1)
    
    ax.set_xlabel("X Axis")
    ax.set_ylabel("Y Axis")
    ax.set_zlabel("Z Axis")
    ax.set_title(str(hyperplane))

    for i in range(5):
        ## seed

        random_start = (random.random() * 2 - 1, random.random() * 2 - 1, random.random() * 2 - 1, random.random() * 2 - 1)
        
        ws, xs, ys, zs, crossings = experiment3(   random_start,
                                                    T = 800, 
                                                    lmbda = [0.086, 0.141, 0.773],
                                                    hyperplane = HyperPlane(1, 1, 1, 1, 1) )

        ## slice
        crossings_array = np.array(crossings)
        indices = list(np.where(crossings_array < 0)[0])
        ws = list(np.array(ws)[indices])
        xs = list(np.array(xs)[indices])
        ys = list(np.array(ys)[indices])
        ax.scatter(ws, xs, ys, label=str(i))

    plt.show()


def P(x_0, hyperplane, lmbda = default_lmbda, dt = 0.01, TMAX = 10000):
    """
    x: a point on C2 intersecting with a hyperplane
    P(x): the first intersection with the hyperplane we experience again
    """
    stepCnt_MAX = math.ceil(TMAX / dt)

    # set initial point
    # x_1 = x_0[0]
    # y_1 = x_0[1]
    # x_2 = x_0[2]
    # y_2 = x_0[3]

    current_pt = list(x_0)
    # track intersections with hyperplane
    intersect_checker = IntersectChecker(hyperplane)
    
    # Stepping through "time".
    for i in range(stepCnt_MAX):
        # Derivatives of the W, X, Y, Z state
        derivs = dots(current_pt, lmbda )

        ## copy last pt
        old_pt = list(current_pt)

        ## compute new point
        for j in range(4):
            current_pt[j] = old_pt[j] + (derivs[j] * dt)

        # print(hyperplane(current_pt) )

        ## return on first intersect
        if intersect_checker(current_pt) != 0:
            ## explain what happened
            
            if DEBUG:
                print("x_0 = " + str(x_0) + ": " + str(hyperplane(x_0)))
                print("P(x_0) = " + str(current_pt) + ": " + str(hyperplane(current_pt)))
            return current_pt

    raise RuntimeError("P was unable to find intersection before step count maximum: " + str(stepCnt_MAX))
    return None

def experiment5():
    """
    Get Poincare intersections, try to minimize P
    """
    hyperplane = HyperPlane(1, 1, 1, 1, 1)
    
    # x_0 = (-0.078439752265828444, -0.11655096717320924, 1.9097262049137886, -0.70963738873898563)
    # print(P(x_0, hyperplane))

    # x_0 = (-0.61168029983645567, -0.12264391918900905, 0.47264337024775721, 1.2702327130662128)
    # print(P(x_0, hyperplane))

    # x_0 = (0.66917299080611725, 1.4403296413593518, -3.012570027015701, 1.9493456489793897)
    # print(P(x_0, hyperplane))

    x_0 = (0.57845893030192974, 0.20956697187490533, 0.16072556774857419, 0.050224047088646813)
    x_0 = default_start
    x_0 = (-1.4273163021928814, 0.11473207108128994, -1.5062853531752913, -0.441733188102787)
    x_0 = (8.959135804139336, -8.088945521150906, 8.733959499529313, -8.302512165691056)
    x_0 = (-0.4478236838274165, 1.4545418485115778, -0.0033497949221810265, -0.005701739181132484)
    print(P(x_0, hyperplane))

import scipy.integrate as integrate
import scipy.special as special
from scipy.integrate import quad
from scipy.optimize import newton

def experiment6(x_0, hyperplane = HyperPlane(1, 1, 1, 1, 1)):
    

    def P_instance(x):
        # print(x)
        # print(P(x_0 = x, hyperplane = hyperplane))
        return quad_sq_distance(x, P(P(x_0 = x, hyperplane = hyperplane), hyperplane = hyperplane) )

    def newton_search(x_0, N = 25):

        x = list(x_0)

        hessian = nd.core.Hessian(P_instance)
        jacobian = nd.core.Jacobian(P_instance)

        print(x)
        print(hessian )

        for i in range(N):    
            print("PASS " + str(i))
            print("x = " + str(x) )
            adjust = np.matmul(np.linalg.inv(hessian(x)), np.transpose( jacobian(x)))
            adjust = np.transpose(adjust)[0]
            
            ##print(adjust)
            print("dist = " + str(P_instance(x)))
            x = list_subtract(x, adjust)

        print("dist = " + str(P_instance(x)))
        print(x)

    
    #x_0 = (0.54595411849672337, -0.51509962495305384, 0.83292233890090173, 0.13541863005974944)

    newton_search(x_0)

def experiment7(x_0, hyperplane = HyperPlane(1, 1, 1, 1, 1)):
    

    def P_instance(x):
        # print(x)
        # print(P(x_0 = x, hyperplane = hyperplane))
        return quad_sq_distance(x, P(x_0 = x, hyperplane = hyperplane) )

    def newton_search(x_0, N = 25):

        x = list(x_0)

        hessian = nd.core.Hessian(P_instance)
        jacobian = nd.core.Jacobian(P_instance)

        print(x)
        print(hessian )

        for i in range(N):    
            print("PASS " + str(i))
            print("x = " + str(x) )
            adjust = np.matmul(np.linalg.inv(hessian(x)), np.transpose( jacobian(x)))
            adjust = np.transpose(adjust)[0]
            
            ##print(adjust)
            print("dist = " + str(P_instance(x)))
            x = list_subtract(x, adjust)

        print("dist = " + str(P_instance(x)))
        print(x)

    
    #x_0 = (0.54595411849672337, -0.51509962495305384, 0.83292233890090173, 0.13541863005974944)

    newton_search(x_0)

if __name__=="__main__":
    experiment3(expmt = 'trace', hyperplane = HyperPlane(1, 1, 1, 1, 1) )
    #experiment5()
    print(P((-1.4273163021928814, 0.11473207108128994, -1.5062853531752913, -0.441733188102787), hyperplane = HyperPlane(1, 1, 1, 1, 1)) )
    
    ## BLOCK 1
    #experiment6(x_0 = (0.54595411849672337, -0.51509962495305384, 0.83292233890090173, 0.13541863005974944) )
    # x_6 = (0.7837303047229669, -0.27442810285562941, 1.0265787472933432, 0.0029171648807551354)
    # x_7 = (0.47146440737188888, -0.27271751341540784, 0.63597539304550921, -0.13796889725689748)
    # x_8 = (18.820041529355215, -12.266087428818587, 0.56948371607381676, -2.0459982765239517)
    # experiment6(x_0 = x_8)

    ## BLOCK 2
    # experiment6(x_0 = default_start)

    # x_9 = (0.1 , 0.1, 0.123, -0.323)
    # experiment6(x_0 = x_9,  hyperplane = HyperPlane(1, 1, 1, 1, 0))

    ### CURRETN
    ## x_9 = (0.1 , 0.1, 0.123, -0.323)
    ## xperiment7(x_0 = x_9,  hyperplane = HyperPlane(1, 1, 1, 1, 0))

    #hyperplane = HyperPlane(1, 1, 1, 1, 1)
    #print( hyperplane((-1.4273163021928814, 0.11473207108128994, -1.5062853531752913, -0.441733188102787)) )

### Interesting start points:
## (-0.984990993377135, -0.639472461280798, 0.9006369573472872, 0.7777707421666811)
## TODO: build results-saver


