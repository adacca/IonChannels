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

source_file = "model=KBZR-MGWMN_prot_name=ap_dt=0.001_comp.jld2" #name of file to be loaded (path not needed, extension needed)
source_dir = "computations/HERG_comp" #source folder in /data
MGWMN = load(datadir(source_dir, source_file))

source_file = "model=KBZR-WLMSR_prot_name=ap_dt=0.001_comp.jld2" #name of file to be loaded (path not needed, extension needed)
source_dir = "computations/HERG_comp" #source folder in /data
WLMSR = load(datadir(source_dir, source_file))

source_file = "model=KBZR-OGD_prot_name=ap_dt=0.001_comp.jld2" #name of file to be loaded (path not needed, extension needed)
source_dir = "computations/HERG_comp" #source folder in /data
OGD = load(datadir(source_dir, source_file))


#user analysis
using Plots
plt1 = plot(legend=:topleft,xaxis="time [ms]",yaxis="Excess heat dissipation [1/kBT]",size=(400,400))
plt1 = plot!(twinx(),CR["tL"][:,1],CR["tL"][:,2],linestyle=:dash,linecolor=:grey,ylabel="Transmembrane voltage [mV]",label=:none)
plt1 = plot!(CR["tL"][:,1],CR["Q_ex"],label="CR")
plt1 = plot!(MGWMN["tL"][:,1],MGWMN["Q_ex"],label="MGWMN")
plt1 = plot!(WLMSR["tL"][:,1],WLMSR["Q_ex"],label="WLMSR")
plt1 = plot!(OGD["tL"][:,1],OGD["Q_ex"],label="OGD")
plt1 = plot!(title="A",titlelocation=:left) #should be size in mm

#display(plt1)
#readline()

#saving data 
figname= "HERG_ap_qex.png"
wsave(datadir("figs",figname),plt1)
