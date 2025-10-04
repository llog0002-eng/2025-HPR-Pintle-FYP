def Dyer_N2O(P1,P2,T1,Cd,A):

    from HEM_N2O_func import HEM_N2O
    from SPI_N2O_func import SPI_N2O
    import numpy as np
    import CoolProp.CoolProp as CP
    
    Pv = CP.PropsSI("P","T",T1,"Q",1,"N2O") # Pa

    ### HEM

    dP = 0.05*1e5
    dPs = np.arange(0.0001,P1-P2+dP,dP)
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

    return mdots_Dyer_choked[-1]
