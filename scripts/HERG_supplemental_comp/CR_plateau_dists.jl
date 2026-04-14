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
source_file = "model=KBZR-CR_prot_name=ap_dt=0.001.jld2" #name of file to be loaded (path not needed, extension needed)
source_dir = "distributions/HERG_dist" #source folder in /data
CR = load(datadir(source_dir, source_file))
#dict keys: ["tL", "Q_ex", "H"]

#saving data 
#figname= "HERG_ap_dH.png"
#wsave(datadir("figs",figname),plt1)

dss = CR["dists"]-CR["steadystates"]

plt1 = plot(legend=:topleft,title="Distance from steady state (CR)",xaxis="time [ms]",yaxis="Occupancy difference [%]",dpi=200)#,size=(400,400),title="B",titlelocation=:left) #full graph
plt1 = plot!(CR["tL"][:,1],dss,xlims=(300,500),labels=["C3" "C2" "C1" "O" "I"],color=:black,linestyles=[:dot :dash :dashdot :dashdotdot :solid])
plt1 = plot!(twinx(),CR["tL"][:,1],CR["tL"][:,2],linestyle=:dash,linecolor=:grey,ylabel="Transmembrane voltage [mV]",label=:none,xlims=(300,500))

#display(plt1)
#readline()

figname= "CR_plateau_dist.png"
wsave(datadir("figs",figname),plt1)
