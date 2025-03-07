#9-state soduim channel from Pal and Gangopadhyay
"""
	Gmatrix(alpha)

calculates transition rate matrix for 9-state Na+ model, uses naparam(alpha)
inputs:
alpha--Float64, transmembrane voltage in mV
returns:
9x9 matrix of transition rates for Na+ channel
""" 
function Gmatrix(alpha)
     a,b = param(alpha) #get param as f(alpha) a,b are 5-element lists
     #print("params ok") #debugging
     u = 1.2 #fitting param
     #P0(A) P1 P2 P3 P4 P5(A) P6(I) P7 P8(I)
     return [-(u^2*a[1]) u^2*a[1]     0            0                   0          0          0          0          0         ;
               b[1]    -(b[1]+u*a[1]) u*a[1]       0                   0          0          0          0          0         ;
               0         u*b[1]     -(u*b[1]+a[1]) a[1]                0          0          0          0          0         ;
               0         0            u^2*b[1]   -(u^2*b[1]+a[2]+a[4]) a[2]       0          a[4]       0          0         ;
               0         0            0            b[2]              -(b[2]+a[3]) a[3]       0          0          0         ;
               0         0            0            0                   b[3]     -(b[3]+a[5]) 0          0          a[5]      ;
               0         0            0            b[4]                0          0        -(b[4]+a[2]) a[2]       0         ;
               0         0            0            0                   0          0          b[2]     -(b[2]+a[3]) a[3]      ;
               0         0            0            0                   0          b[5]       0          b[3]     -(b[5]+b[3]) ]
    end


"""
	param(alpha)

calculates parameters for 9-state Na+ transition rate matrix, called by gna(alpha)
input:
alpha--Float64, transmembrane voltage in mV
returns:
a,b as tuple, each 5-element list *in per ms (ms^-1) (changed from original that wanted it in fucking SECODNVCS

"""
function param(alpha) #returns parameters for na+ model
    a = 1e-3*[4779 5045 1684 19.8 800] #ms^-1, a_0
    b = 1e-3*[10.3 12.1 2360 0 59.8] #ms^-1, b_0, b4 is constraint param, b4(V) = a4*b5/a5
    q = [2.83 3.16 0.077 5.573 0.16] #gating charge 
    d = [0.053 0.5 0.78 0.12 0.33] #fractional elecrtical distance (dim-less)
    kbTe = 24.4 #mV kbT/e with abs temp T
    for i = 1:5
        a[i] = a[i]*exp(q[i]*alpha*d[i]/kbTe) #alpha is V in mV
        b[i] = b[i]*exp(-q[i]*alpha*(1-d[i])/kbTe) #alpha is V in mV
	end
    b[4] = a[4]*b[5]/a[5] #normalizing value for reversibility
    return (a,b)
    end