#5-state potassium channel from Mazhari-greenstein-Winslow-Marban-Nuss via Bett-Zhou-Rasmusson

"""
	Gmatrix(alpha)

calculates transition rate matrix for 5-state K+ model, calls param(alpha)
inputs:
alpha--Float64, transmembrane voltage in mV
returns:
5x5 matrix of transition rates for K+ channel in ms^-1
""" 
function Gmatrix(alpha)
    a1,a2,ai,ai2,b1,b2,bi,phi = param(alpha) #param as f(alpha)
    Kf,Kb = (0.0266,0.1348)
    # C1 C2 C3 O I
    return [-a1   a1     0            0       0   ;
             b1 -(b1+Kf) Kf           0       0   ;
             0    Kb   -(Kb+a2+ai2)   a2      ai2 ;
             0    0      b2         -(b2+ai)  ai  ;
             0    0      phi          bi    -(phi+bi)]
    end
"""
	param(alpha)

calculates parameters for 5-state K+ transition rate matrix, called by gk(alpha)
input:
alpha--Float64, transmembrane voltage in mV
returns:
an,bn--tuple of floats, two parameters needed for gk()
"""
function param(alpha) #returns parameters for k+ model
    a1 = 0.0069*exp(0.0272*alpha)
    b1 = 0.0272*exp(-0.0431*alpha)
    a2 = 0.0218*exp(0.0262*alpha)
    b2 = 0.0009*exp(-0.0269*alpha)
    ai = 0.0622*exp(0.0120*alpha)
    bi = 0.0059*exp(-0.0443*alpha)
    ai2= 1.29e-5*exp(2.71e-6*alpha)
    phi = ai2*bi*b2/(a2*ai)
    return a1,a2,ai,ai2,b1,b2,bi,phi
    end