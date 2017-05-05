# Code to Computer-Assisted Methods for Studying Periodic Orbits on R^4

## Abridged Usage Guide, from Appendix

This is by no means a complete usage guide; we are simply showcasing some principal interfaces here. For additional assistance submit an issue or reach me at yuan@americanmarineresearch.org.

#### IntervalN: 
`IntervalN` implements $n$ dimensional interval vectors. We construct `IntervalN` using one-dimensional intervals provided in PyInterval Taschini. A number of operations are supported. Source located at `src/intervals/intervaln.py`. Usage examples follow:

~~~
# import pyinterval
from interval import *
# import these tools
from intervals.intervaln import intervaln

## construction
x = IntervalN([interval([1, 2]), interval([5, 6]), \
                interval([-2, 4]), interval([-2, -9])])
                
x_clone2 = x.clone() # duplicates 

y = IntervalN([interval([10, 2]), interval([52, 61]), \
                interval([-22, 41]), interval([-34, -10])])

x.__str__() # returns a string representation
str(y) # same as above

## comparisons
x_clone2 == x_clone # returns True
x_clone2 != x_clone # returns False
x_clone2 == y # returns False
x_clone2 != y # returns True

x * 4 # scales componentwise
x in y # return 
x in x * 4 # returns False
x in x * interval([0.999, 1.001]) # returns True

## Operations, behavior as expected
x + y 
x - y
x * -2.1
x / 2
x.hull() ## returns convex hull of interval vector
~~~

#### Evolution: 
The `Evolution` class provides an interface for describing and studying systems in which the user specifies the evolution function. We provide tools for simulation at desired resolutions, as well as visualization tools. Source located at `src/evolution.py`. Usage paired with that of `Experiment`, below.

#### Experiment: 
The `Experiment` class provides an interface for quickly deploying experiments, documenting results, and replicating said experiments. In particular, the class was built to handle the volume of simulations we ran and the challenges in documenting and replication. The `Experiment` class was designed to be used alongside the `Evolution` and `EvolutionInterval` classes. Note that Experiment will create a unique folder and a number of files based on supplied information. Source located at `src/experiment.py`. Usage examples follow:

~~~
from experiment import *
from evolution import *
from thesis_defaults import *

evo = EvolutionValdez(lmbda = lmbda_set_1)
print(evo) # string representation of our system
expmt = Experiment(evo = evo, title = "Valdez System Plot", 
                    descr = "Just simulate and plot")

print(expmt) # string representation of our system

## save plots
plt = evo.gen_plot(default_start, plot_type = 1, 
                    stepCnt = 100000, DEBUG = False)
expmt.savePlot(plt)

plt3 = evo.gen_plot(default_start, plot_type = 0, 
                    stepCnt = 100000, DEBUG = False)
expmt.savePlot(plt3)

## save information to textdump file
expmt.print("Text appended to info file successfully")
~~~

#### EvolutionInterval 
The `EvolutionInterval` class provides an interface for rigorously describing and studying systems in which the user specifies the evolution function. It is the set-valued analog of the `Evolution` class, equipped with interval calculus. Source located at `src/evolutioninterval.py`. Usage examples follow:

~~~
from intervals.intervaln import intervaln
from thesis_defaults import *
from evolutioninterval import *
from experiment import *

## construction
ev = Evolution_Valdez(lmbda = lmbda_set_1)
x_k = listToCol([interval(pt) for pt in default_start])
~~~

We implement the analytically computed Jacobian of $f$ and $\Phi$

~~~
jac_x = ev.nabla(intvl_start)
testnablaPhi = ev.nablaPhi(h = 0.01, x = x_k, p = 2)
~~~

We implement analytically computed $p$-th order derivatives of $f$

~~~
ev.dphidt(x_k, p = 1)
ev.dphidt(x_k, p = 2)
ev.dphidt(x_k, p = 3)
~~~

We also implement analytically computed $p$-th order derivatives of $V$

~~~
ev.dVdt(t = 0, x = x_k, p = 0)
ev.dVdt(t = 0, x = x_k, p = 1)
ev.dVdt(t = 0, x = x_k, p = 2)
~~~

Altogether we can generate our solution iteratively 

~~~
data = self.ev.generate(self.params['start_pt'], 
                        h = None,  
                        p_e = 2,
                        stepCt = self.params['step_ct'] )
~~~


## License 
*We provide the MIT open source license for the tools developed for and presented in this paper.*

Copyright 2017 Yuan Wang

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
