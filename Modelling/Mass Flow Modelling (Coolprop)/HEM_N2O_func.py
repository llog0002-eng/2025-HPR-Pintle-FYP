def HEM_N2O(P1,P2,T1,Cd,A):

    import numpy as np
    import CoolProp.CoolProp as CP

    # From "An investigation into hybrid rocket motors" Masters thesis by Jonas Sømod Ahmed

    h1 = CP.PropsSI("H","P",P1,"Q",0,"N2O") # Assumes saturated liquid
    s1 = CP.PropsSI("S","P",P1,"Q",0,"N2O") # Assumes saturated liquid
    
    s2 = s1 # Isentropic assumption
    rho2 = CP.PropsSI("D","P",P2,"S",s1,"N2O")
    h2 = CP.PropsSI("H","P",P2,"S",s1,"N2O")

    mdot = Cd*rho2*A*np.sqrt(2*(h1-h2))

    return mdot

