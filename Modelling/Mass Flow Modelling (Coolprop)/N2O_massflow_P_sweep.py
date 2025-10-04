from HEM_N2O_func import *
from SPI_N2O_func import *
import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP

P2 = 1e5 # Pa
P1 = 50e5 # Pa
T1 = 290 # K
Cd = 0.7
A = 32e-6 # m2
Pv = CP.PropsSI("P","T",T1,"Q",1,"N2O") # Pa

### HEM

dP = 0.05*1e5
dPs = np.arange(0,P1-P2+dP,dP)
n = len(dPs)
mdots_HEM = np.zeros((n,1))
mdots_HEM_choked = np.zeros((n,1))

for i in range(n):
    P2 = P1 - dPs[i]
    mdots_HEM[i] = HEM_N2O(P1,P2,T1,Cd,A)

    # Enforce choked conditions
    if mdots_HEM[i] < mdots_HEM[i-1]:
        mdots_HEM_choked[i] = mdots_HEM_choked[i-1]
    else:
        mdots_HEM_choked[i] = mdots_HEM[i]

### SPI

mdots_SPI = np.zeros((n,1))

for i in range(n):
    P2 = P1 - dPs[i]
    mdots_SPI[i] = SPI_N2O(P1,P2,T1,Cd,A)

### Dyer

mdots_Dyer = np.zeros((n,1))
mdots_Dyer_choked = np.zeros((n,1))

for i in range(n):
    P2 = P1 - dPs[i]
    k = 1 # Assuming upstream is at saturation pressure
    mdots_Dyer[i] = k/(1+k)*mdots_SPI[i] + 1/(1+k)*mdots_HEM[i]
    if mdots_Dyer[i] < mdots_Dyer[i-1]:
        mdots_Dyer_choked[i] = mdots_Dyer_choked[i-1]
    else:
        mdots_Dyer_choked[i] = mdots_Dyer[i]

plt.plot(dPs,mdots_HEM, linestyle = "dashed",color = "b", label = "HEM")
plt.plot(dPs,mdots_HEM_choked,color="b", label = "HEM, choked")
plt.plot(dPs,mdots_SPI,color = "r", label = "SPI")
plt.plot(dPs,mdots_Dyer, linestyle = "dashed", color = "g", label = "Dyer")
plt.plot(dPs,mdots_Dyer_choked,color = "g", label = "Dyer, choked")

plt.title("Mass flow rate vs Pressure differential, A=%.1f mm2"%(A*1e6) + ", Cd = %.2f"%Cd)
plt.legend()
plt.xlabel("P1-P2 [Pa]")
plt.ylabel("Mass Flow Rate [kg/s]")
plt.show()

print(mdots_Dyer_choked[round((P1-P2)/dP)])