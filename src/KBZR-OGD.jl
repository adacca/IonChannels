#4-state potassium channel from Oehmen-Giles-Demir via Bett-Zhou-Rasmusson **[K+] dependence**

"""
	Gmatrix(alpha)

calculates transition rate matrix for 4-state K+ model, calls param(alpha)
inputs:
alpha--Float64, transmembrane voltage in mV
returns:
4x4 matrix of transition rates for K+ channel in ms^-1
""" 
function Gmatrix(alpha)
    a2,ai,b2,bi = param(alpha) #param as f(alpha)
    Kf,Kb = (0.0176,0.684)
    # C2 C3 O I
    return [-Kf   Kf     0      0 ;
             Kb -(Kb+a2) a2     0 ;
             0    b2   -(b2+ai) ai;
             0    0      bi    -bi]
    end
"""
	param(alpha)

calculates parameters for 4-state K+ transition rate matrix, called by gk(alpha)
**[K+] dependence**
input:
alpha--Float64, transmembrane voltage in mV
returns:
an,bn--tuple of floats, two parameters needed for gk()
"""
function param(alpha) #returns parameters for k+ model
    Kcon = 5.4 #[K+] in mM, can vary, but 5.4 seems to be a central value
    a2 = 0.0787*exp(0.0378(alpha+10))
    b2 = 0.0035*exp(-0.0252(alpha+10))
    ai = 0.264/(Kcon/5.4)^0.4*exp(0.0164(alpha+10))
    bi = 0.0849/(Kcon/5.4)^0.05*exp(-0.0454(alpha+10))
    return a2,ai,b2,bi
    end
