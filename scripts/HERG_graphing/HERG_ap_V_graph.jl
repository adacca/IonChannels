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

#P1=-20mV
#retrieving data
source_file = "model=KBZR-CR_prot_name=ap_dt=0.001_comp.jld2" #name of file to be loaded (path not needed, extension needed)
source_dir = "computations/HERG_comp" #source folder in /data
CR = load(datadir(source_dir, source_file))
#dict keys: ["tL", "Q_ex", "H"]


#user analysis
using Plots
plt1 = plot(CR["tL"][1:200000,1],CR["tL"][1:200000,2],linestyle=:dash,linecolor=:grey,linewidth=2,xlabel="time [ms]",ylabel="Transmembrane voltage [mV]",label=:none,size(400,400),dpi=200)
#plt1 = plot!(title="",titlelocation=:left) #should be size in mm

#display(plt1)
#readline()

#saving data 
figname= "HERG_ap_V.png"
wsave(datadir("figs",figname),plt1)
