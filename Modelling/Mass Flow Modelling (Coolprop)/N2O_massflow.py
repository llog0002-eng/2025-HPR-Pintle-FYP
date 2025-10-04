import numpy as np
import matplotlib.pyplot as plt
from Dyer_N2O_func import Dyer_N2O

P2 = 1e5 # Pa
P1 = 50e5 # Pa
T1 = 290 # K
Cd = 0.7
A = 32e-6 # m2

mdot_Dyer = Dyer_N2O(P1,P2,T1,Cd,A)
print(mdot_Dyer)