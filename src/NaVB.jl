#9-state soduim channel from Vanderberg and Bezanilla
"""
	Gmatrix(alpha)

calculates transition rate matrix for 9-state Na+ model, uses naparam(alpha)
inputs:
alpha--Float64, transmembrane voltage in mV
returns:
9x9 matrix of transition rates for Na+ channel (ms^-1)
""" 
function Gmatrix(alpha)
    y,z,a,b,c,d,f,g,i,j = param(alpha) #get param as f(alpha)
    #C1 C2 C3 C4 C5 I4 I5 I O
    return [-(y)  y    0    0      0    0    0    0    0    ;
              z -(z+y) y    0      0    0    0    0    0    ;
              0   z  -(z+y) y      0    0    0    0    0    ;
              0   0    z  -(z+a+g) a    g    0    0    0    ;
              0   0    0    b    -(b+c) 0    0    0    c    ;
              0   0    0    j      0  -(j+a) a    0    0    ;
              0   0    0    0      0    b  -(b+c) c    0    ;
              0   0    0    0      0    0    d  -(d+i) i    ;
              0   0    0    0      d    0    0    f  -(d+f) ]
end


"""
	param(alpha)

calculates parameters for 9-state Na+ transition rate matrix, called by gna(alpha)
input:
alpha--Float64, transmembrane voltage in mV
returns:
y,z,a,b,c,d,f,g,i,j--tuple, parameters needed for gna()
"""
function param(alpha) #returns parameters for na+ model
    y = 16609*exp(1.50*0.22*alpha/24)*1e-3 #converting from s^-1 to ms^-1
    z = 971*exp(-1.50*0.78*alpha/24)*1e-3
    a = 5750*exp(0.42*0.99*alpha/24)*1e-3
    b = 4325*exp(-0.42*0.01*alpha/24)*1e-3
    c = 15669*exp(1.91*0.75*alpha/24)*1e-3
    d = 1361*exp(-1.91*0.25*alpha/24)*1e-3
    f = 432*exp(0.91*0.001*alpha/24)*1e-3
    i = 4*exp(-0.91*0.999*alpha/24)*1e-3 #yes, order is correct
    g = 770*exp(0.91*0.001*alpha/24)*1e-3
    j = g*i/f
    return y,z,a,b,c,d,f,g,i,j
    end