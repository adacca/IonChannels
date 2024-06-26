using DrWatson
@quickactivate "IonChannels"

include(srcdir("ICTools.jl"))
include(srcdir("Na5.jl"))

dt = 1e-3
trange = [-2 10]
tL = spike(200,dt) #spike
mu0 = steadystate(Gmatrix,tL[1,2])
Nadist = evolvedist(Gmatrix,tL[:,2],dt)

Naqex = Qex(Nadist,Gmatrix,tL[:,2],dt)

#plt1 = fxnprotplot(Nadist,tL[:,2],tL[:,1])
#plt2 = fxnprotplot(Naqex,tL[:,2],tL[:,1])


include(srcdir("K5.jl"))

mu0 = steadystate(Gmatrix,tL[1,2])
Kdist = evolvedist(Gmatrix,tL[:,2],dt)
Kqex= Qex(Kdist,Gmatrix,tL[:,2],dt)

pt3 = fxnprotplot([Naqex Kqex],tL[:,2],tL[:,1])
