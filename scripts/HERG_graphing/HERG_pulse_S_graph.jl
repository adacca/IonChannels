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
source_file = "model=KBZR-CR_prot_name=pulse_dt=0.001_P1=-20_comp.jld2" #name of file to be loaded (path not needed, extension needed)
source_dir = "computations/HERG_comp" #source folder in /data
CR = load(datadir(source_dir, source_file))
#dict keys: ["tL", "Q_ex", "H"]

source_file = "model=KBZR-MGWMN_prot_name=pulse_dt=0.001_P1=-20_comp.jld2" #name of file to be loaded (path not needed, extension needed)
source_dir = "computations/HERG_comp" #source folder in /data
MGWMN = load(datadir(source_dir, source_file))

source_file = "model=KBZR-WLMSR_prot_name=pulse_dt=0.001_P1=-20_comp.jld2" #name of file to be loaded (path not needed, extension needed)
source_dir = "computations/HERG_comp" #source folder in /data
WLMSR = load(datadir(source_dir, source_file))

source_file = "model=KBZR-OGD_prot_name=pulse_dt=0.001_P1=-20_comp.jld2" #name of file to be loaded (path not needed, extension needed)
source_dir = "computations/HERG_comp" #source folder in /data
OGD = load(datadir(source_dir, source_file))

#user analysis
using Plots
plt1 = plot(legend=:topright,size=(400,400),xlims=(0,12),ylims=(0,0.1),dpi=200)
plt1 = plot!(CR["tL"][:,1],CR["S"],label="CR",yaxis=" ",color=:black,linestyle=:dot)
plt1 = plot!(MGWMN["tL"][:,1],MGWMN["S"],label="MGWMN",color=:black,linestyle=:dash)
plt1 = plot!(WLMSR["tL"][:,1],WLMSR["S"],label="WLMSR",color=:black,linestyle=:dashdotdot)
plt1 = plot!(OGD["tL"][:,1],OGD["S"],label="OGD",color=:black,linestyle=:solid)

#P1 = 0mV
source_file = "model=KBZR-CR_prot_name=pulse_dt=0.001_P1=0_comp.jld2" #name of file to be loaded (path not needed, extension needed)
source_dir = "computations/HERG_comp" #source folder in /data
CR = load(datadir(source_dir, source_file))
#dict keys: ["tL", "Q_ex", "H"]

source_file = "model=KBZR-MGWMN_prot_name=pulse_dt=0.001_P1=0_comp.jld2" #name of file to be loaded (path not needed, extension needed)
source_dir = "computations/HERG_comp" #source folder in /data
MGWMN = load(datadir(source_dir, source_file))

source_file = "model=KBZR-WLMSR_prot_name=pulse_dt=0.001_P1=0_comp.jld2" #name of file to be loaded (path not needed, extension needed)
source_dir = "computations/HERG_comp" #source folder in /data
WLMSR = load(datadir(source_dir, source_file))

source_file = "model=KBZR-OGD_prot_name=pulse_dt=0.001_P1=0_comp.jld2" #name of file to be loaded (path not needed, extension needed)
source_dir = "computations/HERG_comp" #source folder in /data
OGD = load(datadir(source_dir, source_file))


plt2 = plot(legend=:none,yaxis="Total entropy production [1/kBT/s]",xlims=(0,12),ylims=(0,0.5),dpi=200)
plt2 = plot!(CR["tL"][:,1],CR["S"],label="CR",color=:black,linestyle=:dot)
plt2 = plot!(MGWMN["tL"][:,1],MGWMN["S"],label="MGWMN",color=:black,linestyle=:dash)
plt2 = plot!(WLMSR["tL"][:,1],WLMSR["S"],label="WLMSR",color=:black,linestyle=:dashdotdot)
plt2 = plot!(OGD["tL"][:,1],OGD["S"],label="OGD",color=:black,linestyle=:solid)

#P1-20mV
source_file = "model=KBZR-CR_prot_name=pulse_dt=0.001_P1=20_comp.jld2" #name of file to be loaded (path not needed, extension needed)
source_dir = "computations/HERG_comp" #source folder in /data
CR = load(datadir(source_dir, source_file))
#dict keys: ["tL", "Q_ex", "H"]

source_file = "model=KBZR-MGWMN_prot_name=pulse_dt=0.001_P1=20_comp.jld2" #name of file to be loaded (path not needed, extension needed)
source_dir = "computations/HERG_comp" #source folder in /data
MGWMN = load(datadir(source_dir, source_file))

source_file = "model=KBZR-WLMSR_prot_name=pulse_dt=0.001_P1=20_comp.jld2" #name of file to be loaded (path not needed, extension needed)
source_dir = "computations/HERG_comp" #source folder in /data
WLMSR = load(datadir(source_dir, source_file))

source_file = "model=KBZR-OGD_prot_name=pulse_dt=0.001_P1=20_comp.jld2" #name of file to be loaded (path not needed, extension needed)
source_dir = "computations/HERG_comp" #source folder in /data
OGD = load(datadir(source_dir, source_file))


plt3 = plot(legend=:none,dpi=200)
plt3 = plot!(CR["tL"][:,1],CR["S"],label="CR",xlims=(0,12),ylims=(0,1.2),color=:black,linestyle=:dot)
plt3 = plot!(MGWMN["tL"][:,1],MGWMN["S"],label="MGWMN",color=:black,linestyle=:dash)
plt3 = plot!(WLMSR["tL"][:,1],WLMSR["S"],label="WLMSR",color=:black,linestyle=:dashdotdot)
plt3 = plot!(OGD["tL"][:,1],OGD["S"],label="OGD",xaxis="time [ms]",yaxis=" ",color=:black,linestyle=:solid)

#combining plots
#plt = plot(plt1,plt2,plt3,layout = (3,1),size=(600,1200))
plt = plot(plt1,plt2,plt3,layout = (3,1),size=(300,600),left_margin=(5,:mm),right_margin=(5,:mm),dpi=200)
plt = plot!(title=["(d)" "(e)" "(f)"],titlelocation=:left) #should be size in mm

#saving data 
figname= "HERG_pulse_S.png"
wsave(datadir("figs",figname),plt)
