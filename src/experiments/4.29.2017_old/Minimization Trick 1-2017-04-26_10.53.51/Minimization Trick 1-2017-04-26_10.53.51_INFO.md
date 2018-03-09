**Timestamp:** 2017-04-26_10.53.51

**Experiment:** Minimization Trick 1

**System:**
x1' = -x1 + x1*x2 + x1*x3 
x2' = 3*x2*x4 - 2*x1*x2 
x3' = x3*x4 - x3*x1 + x2*x1 
x4' = -x3*x4 - 3*x2*x4 + x1 


**Description:** Seeking orbits by minimizing |f_T(x) - x|^2

**Parameters:**

{'T': 1, 'start_pt': (0.032, 0.308, -0.1, -0.5)}

**Dump:**

Running Newton Search, varying t
min_t:
0.500075
F(min_t)
[0.020535158596841026, 0.15700646834929746, -0.076939952913073928, -0.36060167403306459]