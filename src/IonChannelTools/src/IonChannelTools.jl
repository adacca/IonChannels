module IonChannelTools

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

Return steady state distribution(s) for transition rate matrix 'G' with parameter(s) 'alpha'.

'alpha' may be a single value or a list of values.
"""
function steadystate(G,alpha)
    if applicable(start, alpha)
        return steadystate.(G,alpha) #ss for each
	# **type/format might be messed up, original on CoCalc makes it into an array
    else
        return transpose(abs.(nullspace(transpose(G(alpha)))))
	end
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
	fxnprotplot(fxn,protocol,t,multi=False,labels=:none)

Plots one or multiple functions with a protocol function on a a separate y-axis.

If 'multi', then each row of 'fxn' will be plotted (each first-dimension element).

haven't tested properly, might not work
"""
function fxnprotplot(fxn,protocol,t,multi=False,labels=:none)
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
	

end # module IonChannelTools
end