#5-state soduim channel from Mikhael's paper
"""
	Gmatrix(alpha)

calculates transition rate matrix for 5-state Na+ model, uses naparam(alpha)
inputs:
alpha--Float64, transmembrane voltage in mV
returns:
5x5 matrix of transition rates for Na+ channel
""" 
function Gmatrix(alpha)
    am,ah,bm = param(alpha) #get param as f(alpha)
    k1,k2,k3 = 6/25,2/5,3/2 #set constants
    return [-3*am 3*am        0           0        0 ;
             bm -(2*am+bm+k1) 2*am        0        k1;
             0    2*bm      -(am+2*bm+k2) am       k2;
             0    0           3*bm      -(3*bm+k3) k3;
             0    0           ah          0       -ah]
    end


"""
	param(alpha)

calculates parameters for 5-state Na+ transition rate matrix, called by gna(alpha)
input:
alpha--Float64, transmembrane voltage in mV
returns:
a_m,a_h,b_m--tuple, three parameters needed for gna()
"""
function param(alpha) #returns parameters for na+ model
    am = ((alpha+40e-3)/10e-3)/(1-exp(-(alpha+40e-3)/10e-3))
    ah = .07*exp(-(alpha+65e-3)/20e-3)
    bm = 4*exp(-(alpha+65e-3)/18e-3)
    return am,ah,bm
    end