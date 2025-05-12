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

#retrieving data
source_file = "model=Na5_prot_name=pulse_dt=0.001.jld2" #name of file to be loaded (path not needed, extension needed)
source_dir = "distributions" #source folder in /data
dict = load(datadir(source_dir, source_file))
#dict keys: ["tL", "model", "dt", "dists", "gitpatch", "gitcommit", "tbound", "script", "steadystates", "prot_fxn", "prot_name", "Largs"]

#user analysis
using Plots
#=
dist_plt = plot(title="5-state Na+ channel state distribution over time",xaxis="Time [ms]",yaxis="Proportion in state [%]",legend=:topright)
label = ["A3" "A2" "A1" "O" "I"]
for i in 1:5
    dist_plt = plot!(dict["tL"][:,1],dict["dists"][:,i], labels=label[i])
end
dist_plt = plot!(twinx(),dict["tL"][:,1],dict["tL"][:,2], yaxis = "Transmembrane Voltage [mV]", linestyle=:dash,seriescolor=:gray,label=:none)
display(dist_plt)
readline()
figname= replace(source_file, r".jld2$"=>"_dists.png")
wsave(datadir("figs",figname),dist_plt)
=#

#excess heat calculation + plotting
Q_ex = IonChannelTools.Qex(dict["dists"],dict["steadystates"],dict["dt"])

H = []
for i in 1:length(dict["dists"][:,1])
    append!(H, IonChannelTools.H(dict["dists"][i,:]))
end

include(srcdir(dict["model"]*".jl"))
S_sys = IonChannelTools.S_array(IonChannelTools.S_sys,Gmatrix,dict["tL"][:,2],dict["dists"])
S_tot = IonChannelTools.S_array(IonChannelTools.S_tot,Gmatrix,dict["tL"][:,2],dict["dists"])


q_plt = plot(title="5-state Na+ channel Q_ex over time",xaxis="Time [ms]",yaxis="Q_ex [k_B*T]",legend=:topright,xlimits=[0,5],ylimits=[-1,30])
q_plt = plot!(twinx(),dict["tL"][:,1],dict["tL"][:,2], yaxis = "Transmembrane Voltage [mV]", linestyle=:dash,seriescolor=:gray,label=:none,xlimits=[0,5],ylimits=[-100,10])

q_plt = plot!(dict["tL"][:,1],Q_ex,label="Qex",color=:red)
q_plt = plot!(dict["tL"][:,1],IonChannelTools.ddt(Q_ex,dict["dt"]),label="dQex",linestyle=:dash,color=:pink)
#q_plt = plot!(dict["tL"][:,1],IonChannelTools.ddt(H,dict["dt"]),label="dH",linestyle=:dash,color=:lightblue)
#q_plt = plot!(dict["tL"][:,1],IonChannelTools.ddt(H,dict["dt"])+IonChannelTools.ddt(Q_ex,dict["dt"]),label="dH+dQ",color=:blue)
#q_plt = plot!(dict["tL"][:,1],S_sys,label="S_sys",linestyle=:dash,color=:pink)
#q_plt = plot!(dict["tL"][:,1],S_tot,label="S_tot",color=:red)


display(q_plt)
readline()
figname= replace(source_file, r".jld2$"=>"_qex.png")
wsave(datadir("figs",figname),q_plt)


#saving data (OPT INPUT: destination)
#dest_dir = "computations" #name of destination directory in "/data/"
#d = @strdict Q_ex  #prepare these values to save
#merge(dict,d)
#new_name = replace(source_file, r".jld2$"=>"")*"_comp.jld2" #new name for file (removes previous extension)
#@tagsave(datadir(dest_dir, new_name), d; storepatch=true)
#print("Distributions saved to /data/"*dest_dir*"/"*new_name*"\n")