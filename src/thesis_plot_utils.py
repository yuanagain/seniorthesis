## thesis_utils_plot.py
import matplotlib
matplotlib.use('Agg') 
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math

from thesis_utils import *
from thesis_defaults import *
from thesis_poincare_utils import *

def param_string(x_0, T, dt):
    """
    State parmeters for experiment documentation, see evolutions.py
    """
    return "x_0 = " + str(x_0) + "; T = " + str(T) + "; dt = " + str(dt)

def gen_plot(ws, xs, ys, zs, plot_type = 0, txt = "", subtitle_fontsize = 8, top_adjust = 0.80):

    if plot_type == 0:
        print("Plotting Double Plot Quad Viz")
        plt.figure(1)

        plt.subplot(2, 1, 1)
        plt.subplots_adjust(top=top_adjust)
        plt.plot(xs, ws)
        #plt.yscale('linear')
        plt.title('x-y')
        plt.grid(True)
        #plt.gca().set_aspect('equal')

        plt.subplot(2, 1, 2)
        plt.plot(ys, zs)
        #plt.yscale('linear')
        plt.title('w-z')
        plt.grid(True)
        #plt.gca().set_aspect('equal')
        plt.subplots_adjust(top=top_adjust)
        plt.suptitle(txt, fontsize=subtitle_fontsize)

        return plt

    elif plot_type == 1:
        print("Plotting Overlain Double Plot Quad Viz")
        plt.figure(1)

        plt.plot(xs, ws)

        plt.plot(ys, zs)
        #plt.yscale('linear')
        plt.title('x-w, y-z')
        plt.grid(True)
        #plt.gca().set_aspect('equal')
        plt.subplots_adjust(top=top_adjust + 0.02)
        plt.suptitle(txt, fontsize=subtitle_fontsize )

        return plt

    elif plot_type == 2:
        print("Plotting Sphere Plot Quad Viz")
        fig = plt.figure()
        ax = fig.gca(projection='3d')

        plt.subplots_adjust(top=top_adjust)
        plt.suptitle(txt, fontsize=subtitle_fontsize)

        qdist = quad_distance(ws, xs, ys, zs)

        ws = np.divide(ws, qdist)
        xs = np.divide(xs, qdist)
        ys = np.divide(ys, qdist)
        zs = np.divide(zs, qdist)
        
        ax.plot(xs, ys, zs)
        ax.set_xlabel("X Axis")
        ax.set_ylabel("Y Axis")
        ax.set_zlabel("Z Axis")
        ax.set_title("x/d-y/d-z/d")
        

        return plt

    else: 
        print("Invalid Plot Type")

def plot_quad(ws, xs, ys, zs, plot_type = 0, txt = "", subtitle_fontsize = 8, top_adjust = 0.82):  

    plt = gen_plot(ws, xs, ys, zs, plot_type, txt, subtitle_fontsize, top_adjust)
    plt.show()
    return   

   