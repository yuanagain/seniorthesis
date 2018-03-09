**Timestamp:** 2017-04-26_12.20.25

**Experiment:** Minimization Trick 3

**System:**
x1' = -x1 + x1*x2 + x1*x3 
x2' = 3*x2*x4 - 2*x1*x2 
x3' = x3*x4 - x3*x1 + x2*x1 
x4' = -x3*x4 - 3*x2*x4 + x1 


**Description:** Seeking orbits by minimizing |f_T(x) - x|^2

**Parameters:**

{'start_T': 2, 'start_pt': [1, 2, 1, -1]}

**Dump:**

Running Newton Search, Varying x_0
*(x, g(x)):*
([2.5423772616251941, -0.16999957972850499, 0.21664816982969359, 2.2777014466105729, 0.9657777289405991], nan)