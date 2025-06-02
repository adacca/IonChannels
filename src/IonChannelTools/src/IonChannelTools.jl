module IonChannelTools

using LinearAlgebra, Plots, DataFrames

export spike

"""
	evolvedist(G,L,dt,mu0=steadystate(G,L[1]))

Compute state distributions according to the transition rate matrix function 'G'.

# Arguments
- `G`: function which takes a parameter argument (from 'L') and returns a square transition rate matrix
- `L`: list of parameter values to be fed to 'G', (size determines number of steps to take)
- `dt`: time step
- `mu0`: initial distribution as a row vector, uses steady state if unspecified

"""
function evolvedist(G,L,dt,mu0=steadystate(G,L[1]))
    mu = Array{Float64}(undef,0,length(mu0))
    mu = [mu; mu0] #array of state dist. begins with initial mu0
    for i = 1:length(L)
        Ga = G(L[i]) #calculate transition rate matrix for given driving L
        newmu = transpose(mu[end,:]) * exp(Ga*dt) #new state distribution
        c = findall(x -> x < 0, newmu) #indices of all negative values in newmu
        newmu[c] .= 1e-20

        mu = [mu ; newmu/sum(newmu)] #normalizes to sum(prob)=1 so no lil artifacts hopefully
        end
    return mu[2:end,:] #returns all but the first line (so array size matches L)
    end

"""
	quick_evolve(G,t,mu0)

Compute state distribution directly for time-indep. transition rate matrix

# Arguments
- `G`: transition rate matrix
- `t`: time to evolve in ms
- `mu0`: initial state distribution as row vector
"""
function quick_evolve(G,t,mu0)
    return mu0*exp(t*G) #form of solution verified from class notes from freshman ODE course
end

"""
	steadystate(G,alpha)
    steadystate(G::Array)

Calculate steady state distribution for transition rate matrix 'G' with parameter 'alpha'.

'alpha' must be a signle value.

needs to be normalized bc for whatever reason isn't normalized to begin with?? but done now
"""
function steadystate(G,alpha)
    ss = transpose(abs.(LinearAlgebra.nullspace(transpose(G(alpha)))))
    return ss/sum(ss) #normalized
    end
function steadystate(G)
    ss = transpose(abs.(LinearAlgebra.nullspace(transpose(G))))
    return ss/sum(ss) #normalized
    end


"""
	Qex(dists,G,L,dt)
	Qex(dists,ss,dt)

Calculate excess heat for each step in 'dists'.

# Arguments
- `dists`: array of state distributions (as from evolvedist())
- `G`: transition rate matrix function associated with 'dists'
- `L`: list of parameter values 'G' associated with 'dists'
OR
- `ss`: array of hypothetical steady states corresponding to distributions in dists
- `dt`: time step
"""
function Qex(dists::Array,G::Function,L::Array,dt::Float64)
    Qex=[0] #starts at 0
    for i= 1:length(L)-1 #for each step
	newQ = -((transpose(dists[i+1,:]-dists[i,:])*transpose(-log.(steadystate(G,L[i]))))[1])
	#defined equation term
        Qex = [Qex; Qex[end] + newQ] #summation
        end
    return Qex #array of Qex accumulating over the whole interval
    end
function Qex(dists::Array,ss::Array,dt::Float64)
    Qex=[0] #starts at 0
    for i= 1:length(ss[:,1])-1 #for each step
	newQ = -((transpose(dists[i+1,:]-dists[i,:])*(-log.(ss[i,:])))[1])
	#defined equation term
        Qex = [Qex; Qex[end] + newQ] #summation
        end
    return Qex #array of Qex accumulating over the whole interval
    end
"""
	Wex(dists,G,L,dt)
	Wex(dists,ss,dt)

Calculate excess work for each step in 'dists'.

# Arguments
- `dists`: array of state distributions (as from evolvedist())
- `G`: transition rate matrix function associated with 'dists'
- `L`: list of parameter values 'G' associated with 'dists'
OR
- `ss`: array of hypothetical steady states corresponding to distributions in dists
- `dt`: time step
"""
function Wex(dists::Array,G::Function,L::Array,dt::Float64)
    wex=[0] #starts at 0
    for i= 1:length(L)-1 #for each step
        newW = ((transpose(dists[i,:])*transpose(-log.(steadystate(G,L[i+1]))-(-log.(steadystate(G,L[i])))))[1]) #defined equation term
	wex = [wex; wex[end] + newW] #summation
        end
    return wex #array of Wex accumlating over the whole interval
    end
function Wex(dists::Array,ss::Array,dt::Float64)
    wex=[0] #starts at 0
    for i= 1:length(L)-1 #for each step
        newW = ((transpose(dists[i,:])*transpose(-log.(ss[i+1,:])-(-log.([i,:]))))[1]) #defined equation term
	wex = [wex; wex[end] + newW] #summation
        end
    return wex #array of Wex accumlating over the whole interval
    end

"""
    H(mu)

Compute Shannon entropy of state distribution mu.

#Arguments
- `mu`: 1-dim list of single state distribution

NOTE:
computes using natural logarithm. Must match log base used in Q_ex and W_ex calculations
missing a factor of kB/2. doesn't matter really.
"""
function H(mu)
    #shannon entropy of state distribution mu, using nat log
    h = 0
    for x in mu
        h-= x*log(x) #using nat log, could change?
    end
    return h
end

"""
    S_tot(G,dist)

Compute Pal-Gangopadhyay's 'S dot total' total entropy production rate.

#Arguments
- `G`: transition rate matrix of actual values (not the matrix function)
- `mu`: single state distribution 

NOTE: equivalent to sum of S_sys + S_med
"""
function S_tot(G,mu)
    S_t = 0
    len = length(G[1,:]) #side dimension of matrix
    for i in 1:len
        for j in 1:len
            s1 = (G[i,j]*mu[i]-G[j,i]*mu[j])
            s2 = log(abs((G[i,j]*mu[i])/(G[j,i]*mu[j])))
            if isnan(s2)
                s2 = 0
            end
            s = s1*s2
            S_t += s
        end
    end
    return S_t
end

"""
    S_med(G,mu)

Compute Pal-Gangopadhyay's 'S dot medium' entropy flux into environment.

#Arguments
- `G`: transition rate matrix of actual values (not the matrix function)
- `mu`: single state distribution 

"""
function S_med(G,mu)
    S_m = 0
    len = length(G[1,:]) #side dimension of matrix
    for i in 1:len
        for j in 1:len
            s1 = (G[i,j]*mu[i]-G[j,i]*mu[j])
            s2 = log(abs((G[i,j])/(G[j,i])))
            if isnan(s2)
                s2 = 0
            end
            s = s1*s2
            S_m += s
        end
    end
    return S_m
end

"""
    S_sys(G,mu)

Compute Pal-Gangopadhyay's 'S dot system' system entropy production rate.

#Arguments
- `G`: transition rate matrix of actual values (not the matrix function)
- `mu`: single state distribution 
"""
function S_sys(G,mu)
    S_s = 0
    len = length(G[1,:]) #side dimension of matrix
    for i in 1:len
        for j in 1:len
            s1 = (G[i,j]*mu[i]-G[j,i]*mu[j])
            s2 = log(abs((mu[i])/(mu[j])))
            if isnan(s2)
                s2 = 0
            end
            s = s1*s2
            S_s += s
        end
    end
    return S_s
end

"""
    S_array(S_x,Gmatrix,L,dists)

Accumulate epr over time for a list of state distributions

#Arguments
- `S_x`: epr function (S_tot, S_sys, S_med), must take arguments of (transition rate matrix (values), state distribution)
- `Gmatrix`: transition rate matrix function
- `L`: list of parameter values (for G matrix). should match length of dists
- `dists`: array of state distributions

returns epr over time, length of L
"""
function S_array(S_x,Gmatrix,L,dists)
    S = []
    for i in 1:length(L)
        append!(S, S_x(Gmatrix(L[i]),dists[i,:]))
    end
    return S
end

"""
	I_avg(V,V_0,mu,open)

Computes <I> for list of distributions (actually I/open channel conductance but whatever)

# Arguments
- `V`: matrix of driving voltage values in mV corresponding to mu
- `V_0`: membrane rest potential in mV, usually 90mV for Na+ and -80mV for K+
- `mu`: matrix of state distributions
- `open`: row vector indicating open state. ex: [0 0 0 1] for [closed closed closed open]

returns array of I values
"""
function I_avg(V,V_0,mu,open)
    I = Array{Float64}(undef,0,1)
    for i in 1:size(V,1)
        I = [I; (V[i]-V_0)*open*mu[i,:]] #vector mult swapped bc mu default is column
    end
    return I
end

"""
    ddt(x,dt)

Compute first derivative of one-dim list.

#Arguments
- `x`: one-dim list of values
- `dt`: time step between values (or equivalent)

NOTE: returns same size as `x`
"""
function ddt(x,dt)
    dx = Array{Float64}(undef,0,1)
    for i in 1:length(x)
        if i==1
            dx = [dx; (x[i+1]-x[i])/(dt)]
        elseif i==length(x)
            dx = [dx; (x[i]-x[i-1])/(dt)]
        else
            dx = [dx; (x[i+1]-x[i-1])/(2*dt)]
        end
    end
    return dx
end


"""
	fxnprotplot(fxn,protocol,t,multi=false,labels=:none)

Plots one or multiple functions with a protocol function on a a separate y-axis.

If 'multi', then each row of 'fxn' will be plotted (each first-dimension element).

haven't tested properly, might not work
"""
function fxnprotplot(fxn,protocol,t,multi=false,labels=:none)
    plt = plot(xaxis="Time [ms]")
    if multi
	for i in 1:size(fxn,1) #each row
	    plt = plot!(t,fxn[i,:])
	    end
    else
	plt = plot!(t,fxn)
	end
    plt = plot!(twinx(),t,protocol, yaxis = "Transmembrane Voltage [mV]", linestyle=:dash,seriescolor=:gray,label=:none)
    #plots the driving protocol tL on the same graph (but with a seperate y-axis)
    return plt
    end


"""
[time in ms over trange alpha in mV] with pulse param [high start, high end,low,high]
"""
function pulse(trange,dt,args=[0 5 -100 10]) #square pulse
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
	

"""
the spike from Izhikevich
"""
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

    return [t v] #returns in mV
end

end # module IonChannelTools