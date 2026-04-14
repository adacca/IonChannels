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
source_file = "HERG_ap_ssH_comp.jld2" #name of file to be loaded (path not needed, extension needed)
source_dir = "computations/HERG_comp" #source folder in /data
ssH = load(datadir(source_dir, source_file))


using Plots
Hmax = IonChannelTools.H([0.2 0.2 0.2 0.2 0.2])
OGDmax = IonChannelTools.H([0.25 0.25 0.25 0.25])
plt1 = plot(legend=:top,xaxis="time [ms]",yaxis="Relative Shannon Entropy",dpi=200)#,size=(400,400),title="B",titlelocation=:left) #full graph
#plt1 = plot(legend=:topleft,xaxis="time [ms]",yaxis="Total entropy production [1/kBT]",size=(400,400),xlims=(65,80),ylims=(-0.01,0.4))
plt1 = plot!(twinx(),CR["tL"][:,1],CR["tL"][:,2],linestyle=:dash,linecolor=:grey,ylabel="Transmembrane voltage [mV]",label=:none,xlims=(300,500))
plt1 = plot!(CR["tL"][:,1],ssH["crH"]./Hmax,label="CR",xlims=(300,500),color=:black,linestyle=:dot)
plt1 = plot!(CR["tL"][:,1],ssH["mgwmnH"]./Hmax,label="MGWMN",xlims=(300,500),color=:black,linestyle=:dash)
plt1 = plot!(CR["tL"][:,1],ssH["wlmsrH"]./Hmax,label="WLMSR",xlims=(300,500),color=:black,linestyle=:dashdotdot)
plt1 = plot!(CR["tL"][:,1],ssH["ogdH"]./OGDmax,label="OGD",xlims=(300,500),color=:black,linestyle=:solid)
#display(plt1)
#readline()

#saving data 
figname= "HERG_ap_ssH_rel.png"
wsave(datadir("figs",figname),plt1)
