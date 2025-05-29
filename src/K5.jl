#5-state potassium channel from Mikhael's paper

"""
	Gmatrix(alpha)

calculates transition rate matrix for 5-state K+ model, calls kparam(alpha)
inputs:
alpha--Float64, transmembrane voltage in mV
returns:
5x5 matrix of transition rates for K+ channel
""" 
function Gmatrix(alpha)
    an,bn = param(alpha) #param as f(alpha)
    return [-4*an 4*an     0          0        0;
             bn -(3*an+bn) 3*an       0        0;
             0    2*bn   -(2*an+2*bn) 2*an     0;
             0    0        3*bn     -(3*bn+an) an;
             0    0        0          4*bn    -4*bn]
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
    an = ((alpha+55)/100)/(1-exp(-(alpha+55)/10))
    bn = 1/8*exp(-(alpha+65)/80)
    return an,bn
    end
