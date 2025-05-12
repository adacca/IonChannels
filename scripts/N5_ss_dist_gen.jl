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
model = "Na5" #file in "/src/" containing function "Gmatrix" for transition rate matrix. no suffix. path not needed
include(srcdir(model*".jl"))

#INPUT: definition of driving protocol
#dt = 1e-3 #timestep in ms
#tbound = ... #t boundary in ms, max value or [start end] as relevant to driving function
#Largs = ... #arguments for protocol function, array or string as relevant
#prot_fxn = IonChannelTools.... #function for calculating protocol, must take input in form (tbound,dt,args)
prot_name = "ss_range" #name of protocol as str for file labeling
L = range(start= -100.001,step=2e-2, stop=0.02)



#distributions: actual and steady state (INPUT UNCESESSARY)
#dists = IonChannelTools.evolvedist(Gmatrix,tL[:,2],dt) #distributions over protocol
global steadystates = Array{Float64}(undef,0,5)
for i in 1:length(L) #hypothetical steady states over protocol
    global steadystates = [steadystates; IonChannelTools.steadystate(Gmatrix,L[i])]
end

#saving data (OPT INPUT: destination directory)
directory = "distributions" #name of destination directory in "/data/"
d = @strdict model prot_name L steadystates #save these parameters and values as dict
name = savename(d, "jld2"; connector = "_", ignores = ["steadystates"], sort=false)
@tagsave(datadir(directory, name), d; storepatch=true)
print("Distributions saved to /data/"*directory*"/"*name*"\n")
