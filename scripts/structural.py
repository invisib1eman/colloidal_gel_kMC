
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import MDAnalysis as mda
import freud
import shutil
# Define calculation for direct structure factor

def compute_sf_for_timestep_Direct(u,ts_index1,ts_index2,bins,reset=False):
    ts=u.trajectory[ts_index1]
    box = ts.dimensions[:3]
    deltak=4*np.pi/box[0]
    k_max=(bins+0.5)*deltak
    k_min=0.5*deltak
    sfDirect = freud.diffraction.StaticStructureFactorDirect(bins=bins, k_max=k_max, k_min=k_min)
    for ts in range(ts_index1,ts_index2):
        ts = u.trajectory[ts]
        positions = u.atoms.positions.copy()
        sfDirect.compute((box, positions),reset=reset)
    return sfDirect

# Load trajectory


salt_list = [0.9, 1.7, 3, 9, 33, 66, 80]
runstep = 1000
for salt in salt_list:
    
    dataname=f"../data/{salt}mM_Dump.data"
    trjname=f"../trajectories/{salt}mM_Dump.lammpstrj"
    desttrjname=f"../trajectories/{salt}mM_Dump.lammpsdump"
    shutil.copyfile(trjname,desttrjname)
    # read the trajectory
    u = mda.Universe(dataname, desttrjname)
    bins = 30
    start_index = np.arange(0, 1000, 100)
    end_index = np.arange(10, 1000+10, 100)
    for j in range(0,10):
        sfDirect = compute_sf_for_timestep_Direct(u,start_index[j],end_index[j],bins)
        #plot and compare debye and direct
        k_unit = sfDirect.bin_centers / 15
        plt.plot(k_unit, sfDirect.S_k, label=f"{start_index[j]}")
        #plt.plot(sfDebye.k_values, sfDebye.S_k, label="Debye")
        plt.title("Static Structure Factor")
        plt.xlabel("$k (A^{-1})$")
        plt.ylabel("$S(k)$")
        plt.xscale('log')
        plt.yscale('log')
        plt.legend()
        plt.savefig(f"../figures/salt_{salt}_time_{start_index[j]}_end_{end_index[j]}.png")
        with open(f"../structural_factors/salt_{salt}_time_{start_index[j]}_end_{end_index[j]}.txt",'w') as f:
            f.write(f"k_unit    S_k\n")
            for m in range(0,len(k_unit)):
                f.write(f"{k_unit[m]}   {sfDirect.S_k[m]}\n")
    plt.clf()
