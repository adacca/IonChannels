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
#    - protocol: driving function to be called with (tbound,dt,Largs)
# OPTIONAL
# - directory: destination directory in "/data/", default is "/distributions"

#INPUT: model selection
model = "NaVB" #file in "/src/" containing function "Gmatrix" for transition rate matrix. no suffix. path not needed
include(srcdir(model*".jl"))

#INPUT: definition of driving protocol
dt = 1e-3 #timestep in ms
tbound = 200 #t boundary in ms, max value or [start end] as relevant to driving function
Largs = "rs" #arguments for protocol function, array or string as relevant
prot_fxn = IonChannelTools.spike #function for calculating protocol, must take input in form (tbound,dt,args)
prot_name = "spike" #name of protocol as str for file labeling
tL = prot_fxn(tbound,dt,Largs)

#distributions: actual and steady state (INPUT UNCESESSARY)
dists = IonChannelTools.evolvedist(Gmatrix,tL[:,2],dt) #distributions over protocol
global steadystates = Array{Float64}(undef,0,length(dists[1,:]))
for i in 1:size(tL,1) #hypothetical steady states over protocol
    global steadystates = [steadystates; IonChannelTools.steadystate(Gmatrix,tL[i,2])]
end

#saving data (OPT INPUT: destination directory)
directory = "distributions" #name of destination directory in "/data/"
d = @strdict model dt tbound Largs prot_name prot_fxn tL dists steadystates #save these parameters and values as dict
name = savename(d, "jld2"; connector = "_", ignores = ["prot_fxn" "tL" "dists" "steadystates"], sort=false)
@tagsave(datadir(directory, name), d; storepatch=true)
print("Distributions saved to /data/"*directory*"/"*name*"\n")