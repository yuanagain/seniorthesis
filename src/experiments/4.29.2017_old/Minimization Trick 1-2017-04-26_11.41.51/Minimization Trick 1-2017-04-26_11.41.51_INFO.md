**Timestamp:** 2017-04-26_11.41.51

**Experiment:** Minimization Trick 1

**System:**
x1' = -x1 + x1*x2 + x1*x3 
x2' = 3*x2*x4 - 2*x1*x2 
x3' = x3*x4 - x3*x1 + x2*x1 
x4' = -x3*x4 - 3*x2*x4 + x1 


**Description:** Seeking orbits by minimizing |f_T(x) - x|^2

**Parameters:**

{'T': 2, 'start_pt': [-4, 2, -14, 12]}

**Dump:**

Running Newton Search, Varying x_0
*(x, g(x)):*
SUM(X)=-4
([-0.0488644844750592, -1.0673650683935945, 3.1985812112834608, -1.1930040991304978], 5.5745318096005478)
SUM(X)=0.889347559284
([0.049413423009898427, -0.47935469934645702, 1.6417529293789876, 0.71017076925702871], nan)
SUM(X)=1.9219824223
([nan, nan, nan, nan], nan)