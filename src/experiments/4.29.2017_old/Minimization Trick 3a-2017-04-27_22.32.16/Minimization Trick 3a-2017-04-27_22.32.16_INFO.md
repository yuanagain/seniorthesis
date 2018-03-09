**Timestamp:** 2017-04-27_22.32.16

**Experiment:** Minimization Trick 3a

**System:**
Z' = Z * (-e + u_c * C + u_g * G) 
C' = C * (a_c * P - u_c * Z) - m * C * Z 
G' = G * (a_g * P - u_g * Z) + m * C * Z 
P' = P * (-a_g * G - a_c * C) + e * Z 
e = 1.0; u_c = 1.0; u_g = 1.0; a_c = 3.0; a_g = 1.0; m = 1.0

**Description:** Seeking orbits by minimizing |f_T(x) - x|^2

**Parameters:**

{'start_T': 1, 'start_pt': [0.6, 1.0, 1.2, 1.2]}

**Dump:**
