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
source_file = "model=Na5_prot_name=spike_dt=0.001_Largs=rs_tbound=200.jld2" #name of file to be loaded (path not needed, extension needed)
source_dir = "distributions" #source folder in /data
K = load(datadir(source_dir, source_file))
#dict keys: ["tL", "Q_ex", "H"]

#saving data 
#figname= "HERG_ap_dH.png"
#wsave(datadir("figs",figname),plt1)

dss = K["dists"]-K["steadystates"]

plt1 = plot(legend=:topleft,title="Distance from steady state (Na+)",xaxis="time [ms]",yaxis="Occupancy difference [%]")#,size=(400,400),title="B",titlelocation=:left) #full graph
plt1 = plot!(K["tL"][:,1],dss,labels=["C3" "C2" "C1" "O" "I"])
plt1 = plot!(twinx(),K["tL"][:,1],K["tL"][:,2],linestyle=:dash,linecolor=:grey,ylabel="Transmembrane voltage [mV]",label=:none)

display(plt1)
#readline()

figname= "Na_plateau_dist.png"
wsave(datadir("figs",figname),plt1)
