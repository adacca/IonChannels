using LinearAlgebra
using Plots
import Pkg; Pkg.add("LaTeXStrings")
using LaTeXStrings


"""
	evolvestate(mu0,G,L,dt)

Compute state distributions starting at mu0(row vector dim 1xa).
Applies the transition rate matrix(axa) result of G(function of L[i])) with list of parameter values L

y'all make a g that returns an axa matrix, should have a parameter function too that returns a tuple for the params as a function of L[i]
"""
function evolvestate(mu0,G,L,dt)
    mu = Array{Float64}(undef,0,length(mu0))
    mu = [mu; mu0] #array of state dist. begins with initial ss mu0
    for i = 1:length(L)
        Ga = G(L[i]) #calculate transition rate matrix for given driving L
        newmu = transpose(mu[end,:]) * exp(Ga*dt) #new state distribution
        c = findall(x -> x < 0, newmu) #indices of all negative values in newmu
        newmu[c] .= 1e-20

        mu = [mu ; newmu/sum(newmu)] #normalizes to sum(prob)=1 so no lil artifacts hopefully
        end
    return mu[2:end,:] #returns all but the first ss line (so array size matches L)
    end

"""
    evspec(G,a,n=100)

"""
"Calculates eigenvalues of the transition rate matrix G in the range alpha: a=[start end], returns ev (in columns), associated list of alphas"
function evspec(G,a,n=100)
    m = (a[2]-a[1])/n
    ar = lin([0 n],1,[m a[1]]) #this gets a list of alphas to plug in
    ev = non0ev.(G.(ar[:,2])) #gets eigenvalues, ignoring zero
    eva = Array{Float64}(undef,length(ev),length(ev[1])) #pre-allocate array
    for i = 1:length(ev) #reformats for graphing
        eva[i,:]=ev[i]
        end
    return eva,ar[:,2] #eigenvalues and associated alpha values
    end

"gives steady state dist with transitions G with param alpha, opt. parameter range"
function ss(G,alpha,a2=alpha,n=100)
    if alpha==a2 #meaning one alpha, not a range
        return transpose(abs.(nullspace(transpose(G(alpha)))))
    else
        alphas = range(alpha,a2,length=n) #array of alphas
        sss = ss.(G,alphas) #ss for each
        sssa = Array{Float64}(undef,length(alphas),length(sss[1])) #pre-allocates
        for i = 1:length(alphas) #reformatting
            sssa[i,:]=sss[i]
            end
        return sssa
        end
    end

"array:Q_ex for each step in series B with associated time t, matrix function G, and protocol function L (with no t array attached"
function Qex(G,L,dt)
    B = evolvestate(ss(G,L[1]),G,L,dt) #gets the state evolution
    Qex=[0] #starts at 0
    for i= 1:length(L)-1 #for each step
        Qex = [Qex; Qex[end]-((transpose(B[i+1,:]-B[i,:])*transpose(-log.(ss(G,L[i]))))[1])] #defined summation equation
        end
    return Qex #array of Qex accumulating over the whole interval
    end

"array:W_ex for each step in series with matrix function G, and protocol values L"
function Wex(G,L,dt)
    B = evolvestate(ss(G,L[1]),G,L,dt) #gets the state evolution
    wex=[0] #starts at 0
    for i= 1:length(L)-1 #for each step
        wex = [wex; wex[end]+((transpose(B[i,:])*transpose(-log.(ss(G,L[i+1]))-(-log.(ss(G,L[i])))))[1])] #given summation equation
        end
    return wex #array of Wex accumlating over the whole interval
    end



# plotting functions 
"ex=Qex or Wex, needs tuple of Gs, tL is timeProtocol array, returns one plot"
function _explot(ex,Gs,tL,dt,labels=fill(:none,length(Gs)))
    #calculates and plots excess heat or work for each process listed in Gs with driving protocol tL
    Ws = plot()
    for i = 1:length(Gs) #for each model function
        q = ex(Gs[i],tL[:,2],dt) #calculates excess heat or work
        Ws = plot!(tL[:,1],q,label=labels[i]) #adds it to the plot
        end
    Ws = plot!(twinx(),tL[:,1],tL[:,2],linestyle=:dash,seriescolor=:gray,label=:none) #plots the driving protocol tL on the same graph (but with a seperate y-axis)
    return Ws #returns one plot
    end

"returns vector of plots of (log of) each parameter function from Gps(must be input as [], assumes each gives a tuple-like,iterable return), over range alphas, with labels in [[],[],...]"
function paramspecplot(Gps,alphas,labels)
    #this calculates and plots all parameter values over a list of alphas of one or more parameter-returning functions (the ones that feed into the matrix functions)
    all = []
    for i = 1:length(Gps) #loops for each function listed in Gps
        ps = Gps[i].(alphas) #gets a value for each alpha
        psa = Array{Float64}(undef,length(alphas),length(ps[1])) #preallocates
        for j = 1:length(alphas) #reformats
            psa[j,:].=ps[j]
            end
        pl = plot(alphas,log.(psa),label=labels[i]) #makes a plot of the values
        if length(Gps)==1
            return pl
            end
        all = [all; pl] #accumulate a list of plots
        end
    return all #returns list of plots
    end



# supporting functions
"supports evspec: returns non-zero eigenvalues of A"
function non0ev(A)
    b = eigvals(A)
    return deleteat!(b,findall(x-> abs(x)<=1e-12,b))
    end

function con(trange,dt,arg) #also a driving function, returns a constant value
    t = collect(trange[1]:dt:trange[2])
    a = fill(arg,length(t))
    return [t a]
    end
function lin(trange,dt,args) #linear function a = m*t+b
    m,b = args
    t = collect(trange[1]:dt:trange[2])
    a = fill(b,length(t))
    for i = 1:length(t)
        a[i]+=m*t[i]
        end
    return [t a]
    end



# protocol examples
"[time in ms over trange alpha in mV] with pulse param [start,end,low,high]"
function pulse(trange,dt,args=[0 5 -100e-3 10e-3]) #square pulse
    t0,t1,y0,y1=args
    t = collect(trange[1]:dt:trange[2])
    a = fill(y0,length(t))
    for i = 1:length(t)
        if t[i]>t0&&t[i]<t1 #if in the square section, be the other value
            a[i]=y1
            end
        end
    return [t a]
    end

"the spike from Izhikevich"
function spike(tmax,dt,args="rs",)
    #arg types: Regular Spiking, Intrinsically Bursting, CHattering
    C = 100. #transmembrane capacitance in pF

    if args=="m" #mikhael's version from the paper
        vr = -60.
        vt = -40.
        k = 0.7
        a = 0.03
        b = -2.
        c = -50.
        d = 100. #everything in mV and pA
        vpeak = 35. # spike cutoff
        Iin = 80.
        end

    if args=="rs"
        vr = -60.
        vt = -40.
        k = 0.7
        a = 0.03
        b = -2.
        c = -50.
        d = 100. #everything in mV and pA
        vpeak = 35. # spike cutoff
        Iin = 80.
        end

    if args=="ib"
        vr = -75.
        vt = -45.
        k = 1.2
        b = 5.
        vpeak = 50.
        c = -56.
        a = 0.01
        d = 130.
        Iin = 600.
        end

    if args=="ch"
        vr = -75.
        vt = -45.
        k = 1.2
        b = 5.
        vpeak = 50.
        c = -56.
        a = 0.01
        d = 130.
        Iin = 600.
        end

    T = tmax #max time (from 0) in ms
    tau = dt # step in ms
    t = collect(0:tau:T)
    n = length(t) # number of simulation steps

    v = fill(vr, n)
    u = zeros(size(v)) # initial values

    I = zeros(size(v))
    I[1000:end-1000].=Iin # input DC pulse in pA

    for i = 1:n-1 # forward Euler method
        v[i+1] = v[i] + tau*((k*(v[i]-vr)*(v[i]-vt)) + (I[i]-u[i]))/C
        u[i+1] = u[i] + tau*a*(b*(v[i]-vr)-u[i])
        if v[i+1] >= vpeak # a spike is fired!
            v[i] = vpeak # padding the spike amplitude
            v[i+1] = c # membrane voltage reset
            u[i+1] += d # recovery variable update
            end
        end
    #plot(tau*(1:n), v) # plot the result

    return [t v.*1e-3] #returns in mV
    end



# ion channel models
"""
	gna(alpha)

calculates transition rate matrix for 5-state Na+ model, uses naparam(alpha)
inputs:
alpha--Float64, transmembrane voltage in mV
returns:
5x5 matrix of transition rates for Na+ channel
""" 
function gna(alpha)
    am,ah,bm = naparam(alpha) #get param as f(alpha)
    k1,k2,k3 = 6/25,2/5,3/2 #set constants
    return [-3*am 3*am        0           0        0 ;
             bm -(2*am+bm+k1) 2*am        0        k1;
             0    2*bm      -(am+2*bm+k2) am       k2;
             0    0           3*bm      -(3*bm+k3) k3;
             0    0           ah          0       -ah]
    end
"""
	naparam(alpha)

calculates parameters for 5-state Na+ transition rate matrix, called by gna(alpha)
input:
alpha--Float64, transmembrane voltage in mV
returns:
a_m,a_h,b_m--tuple, three parameters needed for gna()
"""
function naparam(alpha) #returns parameters for na+ model
    am = ((alpha+40e-3)/10e-3)/(1-exp(-(alpha+40e-3)/10e-3))
    ah = .07*exp(-(alpha+65e-3)/20e-3)
    bm = 4*exp(-(alpha+65e-3)/18e-3)
    return am,ah,bm
    end

"""
	gna2(alpha)

calculates transition rate matrix for 9-state Na+ model from Pal-Gangopadhyay, calls na2p(alpha)
inputs:
alpha--Float64, transmembrane voltage in mV
returns:
9x9 matrix of transition rates for Na+ channel
"""
function gna2(alpha)
    a,b = na2p(alpha) #each is a list of five parameters f(alpha)
    u = 1.2 #"fitting constant"
    return [-(u^2*a[1]) u^2*a[1] 0 0 0 0 0 0 0; #rows and columns go P0 -> P8
            b[1] -(b[1]+u*a[1]) u*a[1] 0 0 0 0 0 0;
            0 u*b[1] -(u*b[1]+a[1]) a[1] 0 0 0 0 0;
            0 0 u^2*b[1] -(u^2*b[1]+a[2]+a[4]) a[2] 0 a[4] 0 0;
            0 0 0 b[2] -(b[2]+a[3]) a[3] 0 0 0;
            0 0 0 0 b[3] -(b[3]+a[5]) 0 0 a[5];
            0 0 0 b[4] 0 0 -(b[4]+a[2]) a[2] 0;
            0 0 0 0 0 0 b[2] -(b[2]+a[3]) a[3];
            0 0 0 0 0 b[5] 0 b[3] -(b[5]+b[3]);]
    end

"""
	na2p(alpha)

calculates parameters for 9-state Na+ transition rate matrix, called by gna2(alpha)
input:
alpha--Float64, transmembrane voltage in mV
returns:
a,b--tuple of (5,)-vectors, parameters needed for gna2()
"""
function na2p(alpha) #parameters for the P-G na+ transitions
    a0 = [4779 5045 1684 19.8 800] #given
    b0 = [10.3 12.1 2360 0 59.8] #given
    b0[4] = a0[4]*b0[5]/a0[5] #given
    #b0[4] = 19.8*59.8/800
    q = [2.83 3.16 0.077 5.573 0.16] #given
    d = [0.053 0.5 0.78 0.12 0.33] #given
    a = zeros(5) #pre-allocating
    b = zeros(5)
    for i = 1:5 #for each of five alphas and betas
        a[i] = a0[i]*exp(q[i]*alpha*d[i]/24.4e-3)
        b[i] = b0[i]*exp(-q[i]*alpha*(1-d[i])/24.4e-3)
        end
    return a,b
    end

"""
	gk(alpha)

calculates transition rate matrix for 5-state K+ model, calls kparam(alpha)
inputs:
alpha--Float64, transmembrane voltage in mV
returns:
5x5 matrix of transition rates for K+ channel
""" 
function gk(alpha)
    an,bn = kparam(alpha) #param as f(alpha)
    return [-4*an 4*an     0          0        0;
             bn -(3*an+bn) 3*an       0        0;
             0    2*bn   -(2*an+2*bn) 2*an     0;
             0    0        3*bn     -(3*bn+an) an;
             0    0        0          4*bn    -4*bn]
    end
"""
	kparam(alpha)

calculates parameters for 5-state K+ transition rate matrix, called by gk(alpha)
input:
alpha--Float64, transmembrane voltage in mV
returns:
an,bn--tuple of floats, two parameters needed for gk()
"""
function kparam(alpha) #returns parameters for k+ model
    return ((alpha+55e-3)/100e-3)/(1-exp(-(alpha+55e-3)/10e-3)),1/8*exp(-(alpha+65e-3)/80e-3)
    end





#= 
#K square 
tmin,tmax=0,100 #in ms
n = 200000 #number of steps to take (ends up with n+1)
dt = (tmax-tmin)/n

tLpulse = pulse([tmin tmax],dt,[30 70 -70e-3 10e-3])
#GK = tuple(gk)

q = Qex(gk,tLpulse[:,2],dt)

Qpl = plot()
Qpl=plot!(twinx(),tLpulse[:,1],tLpulse[:,2],linestyle=:dash,seriescolor=:gray,label=:none,lw=.8)
Qpl = plot!(tLpulse[:,1],q,seriescolor=:red3,label="K+ channel",lw=1.7)
savefig("QKsquare.svg")

=#

#=
#K spike
tmin,tmax = 0,200
n = 200000
dt = (tmax-tmin)/n

tLspikers = spike(tmax,dt,"rs")
q2 = Qex(gk,tLspikers[:,2],dt)

Qpl2 = plot()
Qpl2=plot!(twinx(),tLspikers[:,1],tLspikers[:,2],linestyle=:dash,seriescolor=:gray,label=:none,lw=.8)
Qpl2 = plot!(tLspikers[:,1],q2,seriescolor=:red3,label="K+ channel",lw=1.7)
savefig("QKspike.svg")
=#


#dist of K over square pulse
tmin,tmax=0,100 #in ms
n = 200000 #number of steps to take (ends up with n+1)
dt = (tmax-tmin)/n

tLpulse = pulse([tmin tmax],dt,[30 70 -70e-3 10e-3])
mu0 = [.2 .2 .2 .2 .2]

dist = evolvestate(mu0,gk,tLpulse[:,2],dt)

distpl = plot(tLpulse[:,1],log.(dist),labels=["A4","A3","A2","A1","O"])
savefig("Ksquare_distribution.svg")