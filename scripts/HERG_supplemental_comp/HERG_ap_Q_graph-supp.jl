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

#want t=300 -> i = 300/dt = 300/1e-3
#need to subtract value at that index from al values
I = Int(300/1e-3)
MGWMN["Q_ex"] .-= MGWMN["Q_ex"][I]
WLMSR["Q_ex"] .-= WLMSR["Q_ex"][I]
OGD["Q_ex"] .-= OGD["Q_ex"][I]

#user analysis
using Plots
plt1 = plot(legend=:topleft,xaxis="time [ms]",yaxis="Heat dissipation [kBT]",title="Total heat dissipation since t=300 ms")#,size=(400,400))
#plt1 = plot!(twinx(),CR["tL"][:,1],CR["tL"][:,2],linestyle=:dash,linecolor=:grey,ylabel="Transmembrane voltage [mV]",label=:none)#,xlims=(0,250))
plt1 = plot!(MGWMN["tL"][I:end,1],MGWMN["Q_ex"][I:end],label="MGWMN",xlims=(300,900))
plt1 = plot!(WLMSR["tL"][I:end,1],WLMSR["Q_ex"][I:end],label="WLMSR",xlims=(300,900))
plt1 = plot!(OGD["tL"][I:end,1],OGD["Q_ex"][I:end],label="OGD",xlims=(300,900))
plt1 = plot!(title="",titlelocation=:left) #should be size in mm
plt1 = plot!(twinx(),MGWMN["tL"][I:end,1],MGWMN["tL"][I:end,2], linestyle=:dash, linecolor=:grey, ylabel="Transmembrane voltage [mV]", label=:none,xlims=(300,900))

#display(plt1)
#readline()

#saving data 
figname= "HERG_qex_shifted.png"
wsave(datadir("figs",figname),plt1)
