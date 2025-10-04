def SPI_N2O(P1,P2,T1,Cd,A):

    import numpy as np
    import CoolProp.CoolProp as CP

    rho = CP.PropsSI("D","P",P1,"T",T1,"N2O") # kg/m3
    mdot_SPI = Cd*A*np.sqrt(2*rho*(P1-P2))
    return mdot_SPI