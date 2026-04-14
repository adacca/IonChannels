using DrWatson
@quickactivate "IonChannels"

using IonChannelTools, DataFrames, JLD2, Plots

#Imports dataframe of distributions from /data/. User can perform desired computations. Saves computed values back to /data/.

#MAKE A COPY OF THIS FILE AND INPUT DATAFRAME AND EXPORT INFORMATION AS NECESSARY

#USER INPUTS
# 
# OPTIONAL
# - source_dir: location in "/data/" of raw dataframe, default is "distributions"
# - dest_dir: destination in "/data/", default 

#P1=-20mV
#retrieving data
source_file = "model=KBZR-WLMSR_prot_name=pulse_dt=0.001_P1=20.jld2" #name of file to be loaded (path not needed, extension needed)
source_dir = "distributions/HERG_dist" #source folder in /data
CR = load(datadir(source_dir, source_file)) #actually wlmsr but idc to change the var name whatever
#dict keys: ["tL", "Q_ex", "H"]

#saving data 
#figname= "HERG_ap_dH.png"
#wsave(datadir("figs",figname),plt1)

dss = CR["dists"]-CR["steadystates"]

plt1 = plot(legend=:topright,title="Distance from steady state (WLMSR)",xaxis="time [ms]",yaxis="Occupancy difference [%]")#,size=(400,400),title="B",titlelocation=:left) #full graph
plt1 = plot!(CR["tL"][:,1],dss,labels=["C3" "C2" "C1" "O" "I"])
#plt1 = plot!(twinx(),CR["tL"][:,1],CR["tL"][:,2],linestyle=:dash,linecolor=:grey,ylabel="Transmembrane voltage [mV]",label=:none)

#display(plt1)
#readline()

figname= "WLMSR_pulse_dist.png"
wsave(datadir("figs",figname),plt1)
