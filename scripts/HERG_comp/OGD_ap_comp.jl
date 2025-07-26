using DrWatson
@quickactivate "IonChannels"

using IonChannelTools, DataFrames, JLD2

#Imports dataframe of distributions from /data/. User can perform desired computations. Saves computed values back to /data/.

#MAKE A COPY OF THIS FILE AND INPUT DATAFRAME AND EXPORT INFORMATION AS NECESSARY

#USER INPUTS
# 
# OPTIONAL
# - source_dir: location in "/data/" of raw dataframe, default is "distributions"
# - dest_dir: destination in "/data/", default 

#retrieving data
source_file = "model=KBZR-OGD_prot_name=ap_dt=0.001.jld2" #name of file to be loaded (path not needed, extension needed)
source_dir = "distributions/HERG_dist" #source folder in /data
dict = load(datadir(source_dir, source_file))
#dict keys: ["tL", "model", "dt", "dists", "gitpatch", "gitcommit", "tbound", "script", "steadystates", "prot_fxn", "prot_name", "Largs"]
include(srcdir(dict["model"]*".jl"))

#user analysis
Q_ex = IonChannelTools.Qex(dict["dists"],dict["steadystates"],dict["dt"])
H = []
for i in length(dict["dists"][:,1])
    append!(H, IonChannelTools.H(dict["dists"][i,:]))
end
I = IonChannelTools.I_avg(dict["tL"][:,2],-80,dict["dists"],[0 0 0 1])

S = IonChannelTools.S_array(IonChannelTools.S_tot,Gmatrix,dict["tL"][:,2],dict["dists"])
tL = dict["tL"]

#saving data (OPT INPUT: destination)
dest_dir = "computations/HERG_comp" #name of destination directory in "/data/"
d = @strdict tL Q_ex I S #prepare these values to save
#merge(dict,d)
new_name = replace(source_file, r".jld2$"=>"_comp.jld2") #new name for file (removes previous extension)
@tagsave(datadir(dest_dir, new_name), d; storepatch=true)
print("Distributions saved to /data/"*dest_dir*"/"*new_name*"\n")