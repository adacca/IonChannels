using DrWatson
@quickactivate "IonChannels"

using IonChannelTools
#include(srcdir("ICTools.jl"))
include(srcdir("Na5.jl"))

dt = 1e-3
trange = [-2 10]
tL = IonChannelTools.spike(200,dt) #spike
mu0 = IonChannelTools.steadystate(Gmatrix,tL[1,2])
Nadist = IonChannelTools.evolvedist(Gmatrix,tL[:,2],dt)

Naqex = IonChannelTools.Qex(Nadist,Gmatrix,tL[:,2],dt)

plt1 = IonChannelTools.fxnprotplot(Nadist,tL[:,2],tL[:,1])
plt2 = IonChannelTools.fxnprotplot(Naqex,tL[:,2],tL[:,1])
display(plt1)
display(plt2)

include(srcdir("K5.jl"))

mu0 = IonChannelTools.steadystate(Gmatrix,tL[1,2])
Kdist = IonChannelTools.evolvedist(Gmatrix,tL[:,2],dt)
Kqex= IonChannelTools.Qex(Kdist,Gmatrix,tL[:,2],dt)

pt3 = IonChannelTools.fxnprotplot([Naqex Kqex],tL[:,2],tL[:,1])
