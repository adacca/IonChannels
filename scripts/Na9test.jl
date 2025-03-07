using DrWatson
@quickactivate "IonChannels"

using IonChannelTools
#include(srcdir("ICTools.jl"))
include(srcdir("Na9.jl"))

dt = 1e-5 #contours exactly with V when 1e-3, trying smaller values
trange = [-2 10] #not used bc we're doing spike rn
#tL = IonChannelTools.spike(200,dt) #spike
#mu0 = IonChannelTools.steadystate(Gmatrix,tL[1,2])
#print(mu0) #debugging

tL = IonChannelTools.pulse(trange,dt)

#println(param(-100e-3))
#println(param(10e-3))
#global as = Array{Float64}(undef,0,5) #setting up empty arrays to fill w params
#global bs = Array{Float64}(undef,0,5)
#for k in tL[:,2] #for each V value
    #newa,newb = param(k) #returns tuple of 5-element matrix
    #global as = [as ; newa]
    #global bs = [bs ; newb]
    #end

#print(as)
#plt3 = IonChannelTools.fxnprotplot(transpose(log.(as)),tL[:,2],tL[:,1],true)
#display(plt3)
#plt4 = IonChannelTools.fxnprotplot(transpose(log.(bs)),tL[:,2],tL[:,1],true)


Nadist = IonChannelTools.evolvedist(Gmatrix,tL[:,2],dt)
Naqex = IonChannelTools.Qex(Nadist,Gmatrix,tL[:,2],dt)

plt1 = IonChannelTools.fxnprotplot(Nadist,tL[:,2],tL[:,1])
display(plt1)
plt2 = IonChannelTools.fxnprotplot(Naqex[1:1000:end],tL[1:1000:end,2],tL[1:1000:end,1])
display(plt2)