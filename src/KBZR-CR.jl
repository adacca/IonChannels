#5-state potassium channel from Clancy-Rudy via Bett-Zhou-Rasmusson **[K+] dependence**

"""
	Gmatrix(alpha)

calculates transition rate matrix for 5-state K+ model, calls param(alpha)
inputs:
alpha--Float64, transmembrane voltage in mV
returns:
5x5 matrix of transition rates for K+ channel in ms^-1
""" 
function Gmatrix(alpha)
    a1,a2,ai,b1,b2,bi,phi = param(alpha) #param as f(alpha)
    Kf,Kb = (2.172,1.077)
    # C1 C2 C3 O I
    return [-a1   a1     0        0      0 ;
             b1 -(b1+Kf) Kf       0      0 ;
             0    Kb   -(Kb+a2+a2)   a2     a2 ;
             0    0      b2     -(b2+ai) ai;
             0    0      phi        bi    -(phi+bi)]
    end
"""
	param(alpha)

calculates parameters for 5-state K+ transition rate matrix, called by gk(alpha)
**[K+] dependence**
input:
alpha--Float64, transmembrane voltage in mV
returns:
an,bn--tuple of floats, two parameters needed for gk()
"""
function param(alpha) #returns parameters for k+ model
    Kcon = 4.5 #[K+] in mM, can vary, but 5.4 seems to be a central value
    a1 = 0.0555*exp(0.05547153(alpha-12))
    b1 = 0.002357*exp(-0.036588*alpha)
    a2 = 0.0655*exp(0.05547153(alpha-36))
    b2 = 0.0029357*exp(-0.02158*alpha)
    ai = 0.656(4.5^0.3/Kcon^0.3)*exp(0.000942*alpha)
    bi = 0.439(4.5/Kcon)*exp(-0.02352(alpha+25))
    phi = bi*b2/ai
    return a1,a2,ai,b1,b2,bi,phi
    end
