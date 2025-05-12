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
source_file = "model=NaVB_prot_name=spike_dt=0.001_Largs=rs_tbound=200.jld2" #name of file to be loaded (path not needed, extension needed)
source_dir = "distributions" #source folder in /data
dict = load(datadir(source_dir, source_file))
#dict keys: ["tL", "model", "dt", "dists", "gitpatch", "gitcommit", "tbound", "script", "steadystates", "prot_fxn", "prot_name", "Largs"]

#user analysis
Q_ex = IonChannelTools.Qex(dict["dists"],dict["steadystates"],dict["dt"])



#saving data (OPT INPUT: destination)
dest_dir = "computations" #name of destination directory in "/data/"
d = @strdict Q_ex  #prepare these values to save
merge(dict,d)
new_name = replace(source_file, r".jld2$"=>"_comp.jld2") #new name for file (removes previous extension)
@tagsave(datadir(dest_dir, new_name), d; storepatch=true)
print("Distributions saved to /data/"*dest_dir*"/"*new_name*"\n")