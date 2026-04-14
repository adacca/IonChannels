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
source_file = "model=Na5_prot_name=spike_dt=0.001_Largs=rs_tbound=200_comp.jld2" #name of file to be loaded (path not needed, extension needed)
source_dir = "computations" #source folder in /data
Na = load(datadir(source_dir, source_file))
#dict keys: ["tL", "Q_ex", "H"]

source_file = "model=K5_prot_name=spike_dt=0.001_Largs=rs_tbound=200_comp.jld2" #name of file to be loaded (path not needed, extension needed)
source_dir = "computations" #source folder in /data
K = load(datadir(source_dir, source_file))


#user analysis
using Plots
plt1 = plot(legend=:topleft,xaxis="time [ms]",yaxis="Heat dissipation [kBT]",dpi=200)#,size=(400,400))
#plt1 = plot!(twinx(),CR["tL"][:,1],CR["tL"][:,2],linestyle=:dash,linecolor=:grey,ylabel="Transmembrane voltage [mV]",label=:none)#,xlims=(0,250))
plt1 = plot!(Na["tL"][:,1],Na["Q_ex"],label="Na+",color=:black,linestyle=:solid)
plt1 = plot!(K["tL"][:,1],K["Q_ex"],label="K+",color=:black,linestyle=:dashdotdot)
plt1 = plot!(title="",titlelocation=:left) #should be size in mm
plt1 = plot!(twinx(),K["tL"][:,1],K["tL"][:,2],linestyle=:dash,linecolor=:grey,ylabel="Transmembrane voltage [mV]",label=:none)

#display(plt1)
#readline()

#saving data 
figname= "nerve_ap_q_rs.png"
wsave(datadir("figs",figname),plt1)
