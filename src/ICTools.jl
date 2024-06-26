using LinearAlgebra, Plots

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
	steadystate(G,alpha)

Return steady state distribution for transition rate matrix 'G' with parameter 'alpha'.

'alpha' must be a signle value.
"""
function steadystate(G,alpha)
    return transpose(abs.(LinearAlgebra.nullspace(transpose(G(alpha)))))
    end


"""
	Qex(dists,G,L,dt)

Calculates excess heat for each step in 'dists'.

# Arguments
- `dists`: array of state distributions (as from evolvedist())
- `G`: transition rate matrix function associated with 'dists'
- `L`: list of parameter values 'G' associated with 'dists'
- `dt`: time step
"""
function Qex(dists,G,L,dt)
    Qex=[0] #starts at 0
    for i= 1:length(L)-1 #for each step
	newQ = -((transpose(dists[i+1,:]-dists[i,:])*transpose(-log.(steadystate(G,L[i]))))[1]) #defined equation term
        Qex = [Qex; Qex[end] + newQ] #summation
        end
    return Qex #array of Qex accumulating over the whole interval
    end

"""
	Wex(dists,G,L,dt)

Calculates excess work for each step in 'dists'.

# Arguments
- `dists`: array of state distributions (as from evolvedist())
- `G`: transition rate matrix function associated with 'dists'
- `L`: list of parameter values 'G' associated with 'dists'
- `dt`: time step
"""
function Wex(dists,G,L,dt)
    wex=[0] #starts at 0
    for i= 1:length(L)-1 #for each step
        newW = ((transpose(dists[i,:])*transpose(-log.(steadystate(G,L[i+1]))-(-log.(steadystate(G,L[i])))))[1]) #defined equation term
	wex = [wex; wex[end] + newW] #summation
        end
    return wex #array of Wex accumlating over the whole interval
    end


"""
	fxnprotplot(fxn,protocol,t,multi=false,labels=:none)

Plots one or multiple functions with a protocol function on a a separate y-axis.

If 'multi', then each row of 'fxn' will be plotted (each first-dimension element).

haven't tested properly, might not work
"""
function fxnprotplot(fxn,protocol,t,multi=false,labels=:none)
    plt = plot()
    if multi
	for i in 1:size(fxn,1) #each row
	    plt = plot!(t,fxn[i,:])
	    end
    else
	plt = plot!(t,fxn)
	end
    plt = plot!(twinx(),t,protocol,linestyle=:dash,seriescolor=:gray,label=:none)
    #plots the driving protocol tL on the same graph (but with a seperate y-axis)
    return plt
    end


"""
[time in ms over trange alpha in mV] with pulse param [start,end,low,high]
"""
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
    ends
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

    return [t v.*1e-3] #returns in mV
    end