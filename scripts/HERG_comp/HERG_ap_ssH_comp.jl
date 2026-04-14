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
source_dir = "distributions/HERG_dist" #source folder in /data

source_file = "model=KBZR-CR_prot_name=ap_dt=0.001.jld2" #name of file to be loaded (path not needed, extension needed)
CR = load(datadir(source_dir, source_file))
#dict keys: ["tL", "Q_ex", "H"]

source_file = "model=KBZR-MGWMN_prot_name=ap_dt=0.001.jld2" #name of file to be loaded (path not needed, extension needed)
MGWMN = load(datadir(source_dir, source_file))

source_file = "model=KBZR-WLMSR_prot_name=ap_dt=0.001.jld2" #name of file to be loaded (path not needed, extension needed)
WLMSR = load(datadir(source_dir, source_file))

source_file = "model=KBZR-OGD_prot_name=ap_dt=0.001.jld2" #name of file to be loaded (path not needed, extension needed)
OGD = load(datadir(source_dir, source_file))

#user analysis
crH = []
for i in 1:length(CR["steadystates"][:,1])
    append!(crH, IonChannelTools.H(CR["steadystates"][i,:]))
end
mgwmnH = []
for i in 1:length(MGWMN["steadystates"][:,1])
    append!(mgwmnH, IonChannelTools.H(MGWMN["steadystates"][i,:]))
end
wlmsrH = []
for i in 1:length(WLMSR["steadystates"][:,1])
    append!(wlmsrH, IonChannelTools.H(WLMSR["steadystates"][i,:]))
end
ogdH = []
for i in 1:length(OGD["steadystates"][:,1])
    append!(ogdH, IonChannelTools.H(OGD["steadystates"][i,:]))
end

#saving data (OPT INPUT: destination)
dest_dir = "computations/HERG_comp" #name of destination directory in "/data/"
d = @strdict crH mgwmnH wlmsrH ogdH #prepare these values to save
#merge(dict,d)
name = "HERG_ap_ssH_comp.jld2"
@tagsave(datadir(dest_dir, name), d; storepatch=true)
print("Distributions saved to /data/"*dest_dir*"/"*name*"\n")