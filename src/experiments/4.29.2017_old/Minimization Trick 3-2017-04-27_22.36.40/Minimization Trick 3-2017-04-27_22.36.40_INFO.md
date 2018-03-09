**Timestamp:** 2017-04-27_22.36.40

**Experiment:** Minimization Trick 3

**System:**
Z' = Z * (-e + u_c * C + u_g * G) 
C' = C * (a_c * P - u_c * Z) - m * C * Z 
G' = G * (a_g * P - u_g * Z) + m * C * Z 
P' = P * (-a_g * G - a_c * C) + e * Z 
e = 1.1; u_c = 1.01; u_g = 1.3; a_c = 3.023; a_g = 1.114; m = 1.121

**Description:** For Colluci Nunez population dynamics system

**Parameters:**

{'start_pt': [0.6, 1.0, 1.2, 1.2], 'start_T': 1}

**Dump:**

Running Newton Search, Varying x_0
*(x, g(x)):*
([0.50382488280726367, 0.83970813801201638, 1.0076497656147423, 1.0076497656147043, 0.76630796428659254], 3.9643046594121012)
([0.40506764813984719, 0.6751127468996867, 0.8101352962798205, 0.81013529627984782, 0.55836842632887751], 1.0749163443896153)
([0.31229569260172657, 0.52049282100281014, 0.62459138520348789, 0.6245913852035071, 0.40061777963275369], 0.29149395348439566)
([0.20823263745368603, 0.34705439575614555, 0.41646527490740803, 0.41646527490740459, 0.29648046444371101], 0.061793424316137463)
([0.15134121464270064, 0.25223535773783812, 0.30268242928540179, 0.30268242928544187, 0.21000598277042176], 0.016154816118987697)
([0.10324870201662259, 0.1720811700277165, 0.20649740403324823, 0.20649740403327718, 0.14535574217994768], 0.0032532452806744222)
([0.076104704645979446, 0.12684117440997478, 0.15220940929196125, 0.15220940929198212, 0.093714553518284449], 0.00074667814964391116)
([0.055077866997306733, 0.0917964449955182, 0.11015573399461569, 0.11015573399462908, 0.060053200925536357], 0.00018591477005258623)
([0.038035025850060251, 0.063391709750106023, 0.076070051700121544, 0.076070051700133701, 0.036610390388263025], 2.8091093123332219e-05)
([0.027176236901250296, 0.045293728168755203, 0.054352473802503687, 0.054352473802511486, 0.022112511828542569], 7.9864893726882309e-06)
([0.022658943546699898, 0.037764905911170765, 0.045317887093402398, 0.045317887093409392, 0.0057078830003213406], 6.0470376223535032e-07)
([0.02698226873144155, 0.044970447885741595, 0.053964537462886361, 0.053964537462894334, 0.0010502409540335744], 8.5747266161072916e-07)
([0.025314199275651886, 0.042190332126094618, 0.050628398551306568, 0.050628398551315373, -0.0030171938177067951], 0.0)