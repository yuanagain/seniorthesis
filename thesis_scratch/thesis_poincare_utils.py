import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math

thesis_poincare_utils_DEBUG = False

class HyperPlane:
    def __init__(self, a, b, c, d, e):
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.e = e
        
    def __call__(self, xywz):
        """Determines which side of the hyperplane that xywz is on"""
        return self.a * xywz[0] + self.b * xywz[1] + self.c * xywz[2] + self.d * xywz[3] - self.e
    
    def whichSide(self, pt):
        val = self.__call__(pt)
        if val > 0: return 1
        elif val < 0: return -1
        else: return 0
        
    def __str__(self):
        return "Hyperplane: " + str(self.a) + "*x_1 + " + \
                        str(self.b) + "*y_1 + " + \
                        str(self.c) + "*x_2 + " + \
                        str(self.d) + "*y_2" + \
                        " = " + str(self.e)
        

class IntersectChecker:
    def __init__(self, hyperplane = HyperPlane(1, 1, 1, 1, 1) ):
        self.hyperplane = hyperplane
        self.flip = 0
        self.last_pt = None
    
    def __call__(self, xywz):
        """
        Checks if we crossed the hyperplane given by abcde. 
        Returns 0 if no crossing. 
        Return -1 if crossed from positive to negative.
        Return 1 if crossed from negative to positive
        """
        val = self.hyperplane.whichSide(xywz)

        

        if self.flip == 0:
            ## oj first pass
            self.flip = val
           
            ## For debugging purposes, check the last pt
            self.last_pt = list(xywz)
            return 0
            
        if val != self.flip:
            ## changed
            if thesis_poincare_utils_DEBUG:
                print("flipped:")
                print(str(self.hyperplane(self.last_pt)) + " -> " + str(self.hyperplane(xywz)) )
            self.flip = val

            ## For debugging purposes, check the last pt
            self.last_pt = list(xywz)
            return val
        
        else:
            ## unchanged

            ## For debugging purposes, check the last pt
            self.last_pt = list(xywz)

            return 0

def poincareExtract(ws, xs, ys, zs, crossings):

    ## slice
    crossings_array = np.array(crossings)
    indices = list(np.where(crossings_array < 0)[0])
    # print("crossings: " + str(len(indices)))
    ws = list(np.array(ws)[indices])
    xs = list(np.array(xs)[indices])
    ys = list(np.array(ys)[indices])
    zs = list(np.array(zs)[indices])
    
    return ws, xs, ys, zs

def poincarePlotCrossings(ws, xs, ys, zs, crossings, txt = " "):
    ## Plot setup
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    plt.subplots_adjust(top=0.85)
    plt.suptitle(txt, fontsize=14)
    
    ## slice
    crossings_array = np.array(crossings)
    indices = list(np.where(crossings_array < 0)[0])
    ws = list(np.array(ws)[indices])
    xs = list(np.array(xs)[indices])
    ys = list(np.array(ys)[indices])
    
    # plt.xlim( (-10, 10) )
    # plt.ylim( (-10, 10) )
    # # plt.zlim( (-10, 10) )
    ax.set_xlim(-3, 3 )
    ax.set_ylim(-3, 3)
    ax.set_zlim(-3, 3)
    ## execute
    ax.scatter(ws, xs, ys)
    ax.set_xlabel("X Axis")
    ax.set_ylabel("Y Axis")
    ax.set_zlabel("Z Axis")
    ax.set_title("x-y-z")

    return plt

def poincarePlot(ws, xs, ys, zs, txt = " ", limits = [None, None, None]):
    ## Plot setup
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    plt.subplots_adjust(top=0.85)
    plt.suptitle(txt, fontsize=14)
    
    # plt.xlim( (-10, 10) )
    # plt.ylim( (-10, 10) )
    # # plt.zlim( (-10, 10) )
    
    if limits[0] != None:
        ax.set_xlim(*limits[0])
    if limits[1] != None:
        ax.set_ylim(*limits[1])
    if limits[2] != None:
        ax.set_zlim(*limits[2])
    # ax.set_ylim(-3, 3)
    # ax.set_zlim(-3, 3)
    ## execute
    ax.scatter(ws, xs, ys)
    ax.set_xlabel("X Axis")
    ax.set_ylabel("Y Axis")
    ax.set_zlabel("Z Axis")
    ax.set_title("x-y-z")

    return plt

def test():
	testplane = HyperPlane(3,2,1,3,-4)
	print(testplane([2,2,2,4,1]))
	print(testplane([-10,-10,-10,-10,-10]))
	print(testplane.whichSide([2,2,2.4,4,1]))
	print(testplane.whichSide([0, 0.0, -4.0, 0, 0]))
	print(testplane.whichSide([-10.3,-10,-10,-10,-10]))
	print(testplane)