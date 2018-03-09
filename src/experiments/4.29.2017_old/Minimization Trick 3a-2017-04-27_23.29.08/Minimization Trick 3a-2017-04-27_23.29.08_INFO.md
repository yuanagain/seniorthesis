**Timestamp:** 2017-04-27_23.29.08

**Experiment:** Minimization Trick 3a

**System:**
Z' = Z * (-e + u_c * C + u_g * G) 
C' = C * (a_c * P - u_c * Z) - m * C * Z 
G' = G * (a_g * P - u_g * Z) + m * C * Z 
P' = P * (-a_g * G - a_c * C) + e * Z 
e = 1.0; u_c = 1.0; u_g = 1.0; a_c = 3.0; a_g = 1.0; m = 1.0

**Description:** For Colluci Nunez population dynamics system

**Parameters:**

{'start_pt': [1, 12, 1, 6], 'start_T': 4}

**Dump:**

Running Newton Search, Varying x_0
*(x, g(x)):*
([3.8822544673241794, 4.831411569987262, 1.9329430304875075, 15.656447081053159, 3.7335534551893921], 64.336588714041255)
([10.54165275135362, -0.97198660937764725, 0.54694001629161981, 14.211065615663205, 3.5819861797564596], nan)