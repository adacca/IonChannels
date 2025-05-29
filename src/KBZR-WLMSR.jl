#5-state potassium channel from Wang-Lu-Morales-Strauss-Rasmusson via Bett-Zhou-Rasmusson

"""
	Gmatrix(alpha)

calculates transition rate matrix for 5-state K+ model, calls param(alpha)
inputs:
alpha--Float64, transmembrane voltage in mV
returns:
5x5 matrix of transition rates for K+ channel
""" 
function Gmatrix(alpha)
    a1,a2,ai,b1,b2,bi = param(alpha) #param as f(alpha)
    Kf,Kb = (0.023761,0.036778)
    # C1 C2 C3 O I
    return [-a1   a1     0        0      0 ;
             b1 -(b1+Kf) Kf       0      0 ;
             0    Kb   -(Kb+a2)   a2     0 ;
             0    0      b2     -(b2+ai) ai;
             0    0      0        bi    -bi]
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
    a1 = 0.022348*exp(0.01176*alpha)
    b1 = 0.047002*exp(-0.0631*alpha)
    a2 = 0.013733*exp(0.038198*alpha)
    b2 = 0.0000689*exp(-0.04178*alpha)
    ai = 0.090821*exp(0.023391*alpha)
    bi = 0.006497*exp(-0.03268*alpha)
    return a1,a2,ai,b1,b2,bi
    end
