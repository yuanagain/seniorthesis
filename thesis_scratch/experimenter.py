## Experiment Runner
## This module helps us track the parameters used in various experiemnts

import os
import sys
import datetime

datetime.datetime.now().strftime("%Y-%m-%d__%H-%M-%S")

class Experiment:
	def __init__(self, lmbda, start_pt, hyperplane = None):
		self.lmbda = lmbda
		self.start_pt = start_pt
		self.hyperplane = hyperplane


	def __str__(self):
		"""returns string representation of self"""
		return "not yet implemented"

	def plot(self):
		"""returns the plots"""
		return

	def run(self):
		"""Run and store results"""
		return

	def export_headers(self):
		""" Save parameters to file"""

		return

def batchExperiments(params, Experiment):

