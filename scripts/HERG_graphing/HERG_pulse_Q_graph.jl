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
plt1 = plot(legend=:topleft,yaxis="h",size=(400,400))
plt1 = plot!(CR["tL"][:,1],CR["Q_ex"],label="CR",yaxis=" ")
plt1 = plot!(MGWMN["tL"][:,1],MGWMN["Q_ex"],label="MGWMN")
plt1 = plot!(WLMSR["tL"][:,1],WLMSR["Q_ex"],label="WLMSR")
plt1 = plot!(OGD["tL"][:,1],OGD["Q_ex"],label="OGD")

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


plt2 = plot(legend=:none,yaxis="Excess heat dissipation [1/kBT]")
plt2 = plot!(CR["tL"][:,1],CR["Q_ex"],label="CR")
plt2 = plot!(MGWMN["tL"][:,1],MGWMN["Q_ex"],label="MGWMN")
plt2 = plot!(WLMSR["tL"][:,1],WLMSR["Q_ex"],label="WLMSR")
plt2 = plot!(OGD["tL"][:,1],OGD["Q_ex"],label="OGD")

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


plt3 = plot(legend=:none)
plt3 = plot!(CR["tL"][:,1],CR["Q_ex"],label="CR")
plt3 = plot!(MGWMN["tL"][:,1],MGWMN["Q_ex"],label="MGWMN")
plt3 = plot!(WLMSR["tL"][:,1],WLMSR["Q_ex"],label="WLMSR")
plt3 = plot!(OGD["tL"][:,1],OGD["Q_ex"],label="OGD",xaxis="time [ms]",yaxis=" ")

#combining plots
#plt = plot(plt1,plt2,plt3,layout = (3,1),size=(600,1200))
plt = plot(plt1,plt2,plt3,layout = (3,1),size=(300,600),left_margin=(5,:mm),right_margin=(5,:mm))
plt = plot!(title=["A" "B" "C"],titlelocation=:left) #should be size in mm

#saving data 
figname= "HERG_pulse_qex.png"
wsave(datadir("figs",figname),plt)
