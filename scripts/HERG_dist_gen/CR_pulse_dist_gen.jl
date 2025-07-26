using DrWatson
@quickactivate "IonChannels"

using IonChannelTools

#Calculates actual and hypothetical steady state distributions of a model with a driving protocol
#and saves model, protocol, and distribution information to a dataframe in /data

#MAKE A COPY OF THIS FILE AND INPUT MODEL AND PROTOCOL INFORMATION AS NECESSARY

#USER INPUTS
# - model: name of file with model
# - protocol definitions:
#    - dt: timestep in ms
#    - tbound: max or range of time domain
#    - Largs: additional arguments for protocol function
#    - prot_fxn: driving function to be called with (tbound,dt,Largs)
#    - prot_name: string name for protocal function for file labeling
# OPTIONAL
# - directory: destination directory in "/data/", default is "/distributions"

#INPUT: model selection
model = "KBZR-CR" #file in "/src/" containing function "Gmatrix" for transition rate matrix. no suffix. path not needed
include(srcdir(model*".jl"))

#INPUT: definition of driving protocol
V0 = -90 #initial holding voltage (in mV)
P1 = 0 #pulse voltage
#P1 = -20 #pulse voltage
#P1 = 20 #pulse voltage
P2 = -40 #voltage after pulse
t1 = 5*1e3 #time at P1 voltage in ms (5s)
t2 = 200 #time at p2 voltage
dt = 1e-3 #timestep in ms
prot_name = "pulse" #name of protocol as str for file labeling
tL = IonChannelTools.pulse([0 t2],dt,[0 0 P2 0])
mu0 = IonChannelTools.steadystate(Gmatrix,V0) #initial ss
mu1 = IonChannelTools.quick_evolve(Gmatrix(P1),t1,mu0) #dist after 5s pulse

#distributions: actual and steady state (INPUT UNCEdSESSARY)
dists = IonChannelTools.evolvedist(Gmatrix,tL[:,2],dt,mu1) #distributions over protocol
global steadystates = Array{Float64}(undef,0,length(dists[1,:]))
for i in 1:size(tL,1) #hypothetical steady states over protocol
    global steadystates = [steadystates; IonChannelTools.steadystate(Gmatrix,tL[i,2])]
end

#saving data (OPT INPUT: destination directory)
directory = "distributions/HERG_dist" #name of destination directory in "/data/"
d = @strdict model dt P1 prot_name tL dists steadystates #save these parameters and values as dict
name = savename(d, "jld2"; connector = "_", ignores = ["tL" "dists" "steadystates"], sort=false)
@tagsave(datadir(directory, name), d; storepatch=true)
print("Distributions saved to /data/"*directory*"/"*name*"\n")