
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import MDAnalysis as mda
import freud
import shutil
import argparse
import os

# Define calculation for direct structure factor
def compute_sf_for_timestep_Direct(u: mda.Universe, ts_index1: int, ts_index2: int, bins: int, reset: bool = False) -> freud.diffraction.StaticStructureFactorDirect:
    """
    Compute structure factor for a given time window using direct method
    
    Args:
        u: MDAnalysis Universe object
        ts_index1: Start timestep index
        ts_index2: End timestep index
        bins: Number of bins for structure factor calculation
        reset: Whether to reset calculation between frames
        
    Returns:
        freud.diffraction.StaticStructureFactorDirect object
    """
    ts = u.trajectory[ts_index1]
    box = ts.dimensions[:3]
    # make deltak finer: np.pi/box[0]
    deltak = np.pi/box[0]
    k_min = 4*deltak
    k_max = (bins)*deltak + k_min
    
    
    sfDirect = freud.diffraction.StaticStructureFactorDirect(bins=bins, k_max=k_max, k_min=k_min)
    
    for ts in range(ts_index1, ts_index2):
        ts = u.trajectory[ts]
        positions = u.atoms.positions.copy()
        sfDirect.compute((box, positions), reset=reset)
    
    return sfDirect
def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Calculate structure factor for a salt concentration')
    parser.add_argument('--salt', type=float, required=True,
                        help='Salt concentration in mM, e.g. 0.9')
    parser.add_argument('--step_size', type=float, required=True,
                       help='Step size for simulation')
    parser.add_argument('--n_particles', type=int, required=True, 
                       help='Number of particles')
    parser.add_argument('--box_length', type=float, required=True,
                       help='Box length L')
    parser.add_argument('--total_time', type=float, required=True,
                       help='Total simulation time')
    parser.add_argument('--tau_crawl', type=float, required=True,
                       help='Tau crawl parameter')
    parser.add_argument('--freeroll', type=bool, required=True,
                       help='Whether the simulation is freeroll')
    
    return parser.parse_args()
def load_trajectory(salt: float, step_size: float, n_particles: int, box_length: float, total_time: float, tau_crawl: float, freeroll: bool) -> mda.Universe:
    """Load trajectory files and return MDAnalysis Universe"""
    if salt == int(salt):
        salt = int(salt)
    if freeroll:
        dataname = f"../../data/{salt}mM_freeroll_step{step_size}_N_{n_particles}_L_{box_length:.0f}_t_{total_time:.0f}_taucrawl_{tau_crawl:.0f}.data"
        trjname = f"../../trajectories/{salt}mM_freeroll_step{step_size}_N_{n_particles}_L_{box_length:.0f}_t_{total_time:.0f}_taucrawl_{tau_crawl:.0f}.lammpstrj"
        desttrjname = f"../../trajectories/{salt}mM_freeroll_step{step_size}_N_{n_particles}_L_{box_length:.0f}_t_{total_time:.0f}_taucrawl_{tau_crawl:.0f}.lammpsdump"
        if os.path.exists(trjname):
            shutil.move(trjname, desttrjname)
    return mda.Universe(dataname, desttrjname)

def save_structure_factors(u: mda.Universe, salt: float, step_size: float, n_particles: int, box_length: float, total_time: float, tau_crawl: float, freeroll: bool, bins: int = 100) -> None:
    """Calculate, plot and save structure factors for different time windows"""
    start_index = np.arange(0, 1000, 100)
    end_index = np.arange(10, 1000+10, 100)
    if salt == int(salt):
        salt = int(salt)
    if freeroll:
        nameprefix = f"{salt}mM_freeroll_step{step_size}_N_{n_particles}_L_{box_length:.0f}_t_{total_time:.0f}_taucrawl_{tau_crawl:.0f}"
    # Create DataFrame to store structure factors
    sf_data = []
    
    for j in range(0, 10):
        sfDirect = compute_sf_for_timestep_Direct(u, start_index[j], end_index[j], bins)
        k_unit = sfDirect.bin_centers / 15
        
        # Store data
        sf_data.append({
            'timestep': start_index[j],
            'k_values': k_unit,
            'S_k': sfDirect.S_k
        })
        
        plt.plot(k_unit, sfDirect.S_k, label=f"t={start_index[j]}")
        
        plt.title(f"Static Structure Factor - {salt}mM")
        plt.xlabel("$k (A^{-1})$")
        plt.ylabel("$S(k)$")
        plt.xscale('log')
        plt.yscale('log')
        plt.legend()
    
        # Save figure
        plt.savefig("../../figures/"+nameprefix+f"_start{start_index[j]}_end{end_index[j]}_finer.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Save structure factor data
    for i, data in enumerate(sf_data):
        df = pd.DataFrame({
            'k_values': data['k_values'],
            'S_k': data['S_k']
        })
        df.to_csv("../../structural_factors/"+nameprefix+f"_start{start_index[i]}_end{end_index[i]}_finer.csv", index=False)



def main():
    """Main function to run structure factor analysis"""
    args = parse_arguments()
    
    # Create output directories if they don't exist
    for dir_path in ['../../figures', '../../data']:
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
    
    # Load trajectory
    u = load_trajectory(args.salt, args.step_size, args.n_particles, args.box_length, args.total_time, args.tau_crawl, args.freeroll)
    
    # Calculate, plot and save structure factors
    save_structure_factors(u, args.salt, args.step_size, args.n_particles, args.box_length, args.total_time, args.tau_crawl, args.freeroll)

if __name__ == "__main__":
    main()



