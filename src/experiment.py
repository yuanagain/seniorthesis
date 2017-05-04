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

import os
import datetime

class Experiment:
    """
    Simply plot
    """
    def __init__(self, evo, title = "TITLE", descr = ""):
        self.title = title
        self.descr = descr
        self.evo = evo
        self.setup()
        self.fig_ct = 0
        self.params = dict()

    def saveParams(self):
        self.saveText("\n**Parameters:**\n")
        self.saveText(str(self.params))
        self.saveText("\n**Dump:**\n")

    def header(self): 
        return "**Experiment:** " + self.title + "\n\n**System:**\n" + str(self.evo) + '\n\n**Description:** ' + self.descr

    def __str__(self):
        return "Experiment: " + self.title + "\n\nSystem:\n" + str(self.evo) + '\n\nDescription: ' + self.descr

    def setup(self):
        """
        Create folder for experiment
        """
        ## experiment and date
        date_txt = datetime.datetime.now().strftime("%Y-%m-%d_%H.%M.%S")
        folder_name = date_txt + "_" + self.title

        pathname = './experiments/' + folder_name
        if not os.path.exists(pathname):
            os.makedirs(pathname)

        self.folder_name = folder_name
        self.date_txt = date_txt
        self.pathname = pathname
        self.infofilename = './' + pathname + '/' + self.folder_name +'_INFO.md'

        
        fle = open(self.infofilename, 'w')

        fle.write("**Timestamp:** " + date_txt + '\n\n')
        fle.write(self.header() )
        fle.close()


        print("Experiment " + self.title + " recorded at " + pathname)

    def savePlot(self, plt):
        """
        Save a plot as png to the directory
        """
        ## increment
        self.fig_ct = self.fig_ct + 1
        plt.savefig(self.pathname + '/' + self.folder_name + '_FIG' + str(self.fig_ct) + '.png')

    def saveText(self, txt):
        """
        Add text to the info folder.
        """
        with open(self.infofilename, 'a') as fle:

            fle.write("\n")
            fle.write(txt)
            fle.close()

    def print(self, txt):
        self.saveText(str(txt))
        print(str(txt))


def main():
    """
    Testing
    """

    print("============")
    evo = Evolution_1a(lmbda = lmbda_set_1)
    print(evo)
    print("============")
    print("============")
    expmt = Experiment(evo = evo, title = "Test-title", descr = "test description")
    print("============")
    print(expmt)

    plt = evo.gen_plot(default_start, plot_type = 0)
    expmt.savePlot(plt)

    plt2 = evo.gen_plot(default_start, plot_type = 1)
    expmt.savePlot(plt2)

    plt3 = evo.gen_plot(default_start, plot_type = 2)
    expmt.savePlot(plt3)

    expmt.saveText("Text appended to info file successfully")

if __name__=="__main__":
    main()
    
