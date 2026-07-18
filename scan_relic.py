import os
import sys
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar

# Path to the compiled C wrapper executable inside the micrOMEGAs model directory
EXECUTABLE = "./calc_relic"

def get_relic_density(mchi, mZp, gp):
    """
    Calls the compiled micrOMEGAs C wrapper to calculate relic density.
    Returns (omega, sigmav_today).
    """
    if not os.path.exists(EXECUTABLE):
        raise FileNotFoundError(f"Executable {EXECUTABLE} not found. Please compile it first.")

    try:
        # Run C wrapper and capture output
        result = subprocess.run(
            [EXECUTABLE, f"{mchi:.8f}", f"{mZp:.8f}", f"{gp:.8e}"],
            capture_output=True,
            text=True,
            check=True
        )
        
        omega = None
        sigmav_today = None
        for line in result.stdout.split('\n'):
            if line.startswith("OMEGA_H2="):
                omega = float(line.split("=")[1])
            elif line.startswith("SIGMAV_TODAY="):
                sigmav_today = float(line.split("=")[1])
                
        if omega is None:
            raise ValueError("Could not parse OMEGA_H2 from output.")
            
        return omega, sigmav_today
    except Exception as e:
        # Return high value on failure (signifies over-abundance / out of bounds)
        return 1e10, 0.0

def objective_function(log_gp, mchi, mZp):
    """
    Objective function for Brent's method: log10(Omega_h2) - log10(0.12)
    We solve in log10(gp) space for better numerical stability.
    """
    gp = 10**log_gp
    omega, _ = get_relic_density(mchi, mZp, gp)
    return np.log10(omega) - np.log10(0.12)

def find_gp_for_relic(mchi, mZp):
    """
    Finds the coupling gp that satisfies Omega h^2 = 0.12 using Brent's method.
    We search for log10(gp) in the interval [-8.0, 0.3] (gp in [1e-8, 2.0]).
    """
    try:
        # Verify brackets have opposite signs
        f_min = objective_function(-8.0, mchi, mZp)
        f_max = objective_function(0.3, mchi, mZp)
        
        if f_min * f_max > 0:
            # If both are positive, even gp=2.0 is over-abundant (no solution)
            if f_max > 0:
                return None
            # If both are negative, even gp=1e-8 is under-abundant
            else:
                return 1e-8
                
        sol = root_scalar(
            objective_function, 
            args=(mchi, mZp), 
            bracket=[-8.0, 0.3], 
            method='brentq', 
            xtol=1e-4
        )
        return 10**sol.root
    except Exception as e:
        print(f"Error finding root for mchi={mchi}, mZp={mZp}: {e}")
        return None

def run_adaptive_scan(mchi, mZp_min=0.001, mZp_max=100000.0, tol=0.15, max_depth=6):
    """
    Performs an adaptive grid scan of mZp to resolve the resonance peak (mZp ~ 2*mchi).
    Starts with a coarse logarithmic grid, then refines regions where gp changes rapidly.
    """
    print(f"\n--- Starting adaptive scan for mchi = {mchi} GeV ---")
    
    # 1. Initialize coarse logarithmic grid (45 points to account for wider scale)
    mZp_grid = np.logspace(np.log10(mZp_min), np.log10(mZp_max), 45)
    
    # Evaluate initial points
    results = []
    for mZp in mZp_grid:
        gp = find_gp_for_relic(mchi, mZp)
        if gp is not None:
            results.append((mZp, gp))
            print(f"mZp: {mZp:.4f} GeV | gp: {gp:.2e}")
            
    # 2. Refine grid iteratively
    depth = 0
    while depth < max_depth:
        print(f"Refinement depth {depth + 1}/{max_depth}...")
        refined = False
        new_results = [results[0]]
        
        for i in range(len(results) - 1):
            p1 = results[i]
            p2 = results[i+1]
            
            log_gp_ratio = abs(np.log10(p1[1]) - np.log10(p2[1]))
            log_mZp_diff = abs(np.log10(p1[0]) - np.log10(p2[0]))
            
            # If gp changes significantly, and the spacing isn't too small, insert midpoint
            if log_gp_ratio > tol and log_mZp_diff > 0.01:
                mZp_mid = 10**((np.log10(p1[0]) + np.log10(p2[0])) / 2.0)
                gp_mid = find_gp_for_relic(mchi, mZp_mid)
                
                if gp_mid is not None:
                    new_results.append((mZp_mid, gp_mid))
                    refined = True
                    print(f"  [Refined] mZp: {mZp_mid:.4f} GeV | gp: {gp_mid:.2e}")
                    
            new_results.append(p2)
            
        results = sorted(new_results, key=lambda x: x[0])
        if not refined:
            print("No further refinement needed.")
            break
        depth += 1
        
    return np.array(results)

def main():
    # DM masses to scan (from Figure 2 of the paper)
    mchi_targets = [0.05, 5.0, 500.0]
    colors = {0.05: 'red', 5.0: 'deepskyblue', 500.0: 'darkviolet'}
    
    plt.figure(figsize=(8, 6))
    
    for mchi in mchi_targets:
        data = run_adaptive_scan(mchi)
        
        # Save raw data to text files
        filename = f"scan_data_mchi_{mchi}.txt"
        np.savetxt(filename, data, header="mZp(GeV) gp", fmt="%.6e")
        print(f"Saved results to {filename}")
        
        # Plot curves
        plt.plot(data[:, 0], data[:, 1], label=f"$m_\\chi = {mchi}$ GeV", color=colors[mchi], lw=2)

    # Styling plot to match Figure 2
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(0.001, 100000.0)
    plt.ylim(1e-8, 1.0)
    plt.xlabel("$m_{Z'}$ (GeV)", fontsize=12)
    plt.ylabel("$g'$", fontsize=12)
    plt.title("$U(1)_{L_\\mu - L_\\tau}$ Relic Abundance Constraint ($\Omega_\\chi h^2 = 0.12$)", fontsize=12)
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend(loc="lower right", fontsize=11)
    
    plt.tight_layout()
    plt.savefig("relic_density_scan.png", dpi=300)
    print("Plot saved to relic_density_scan.png")
    plt.show()

if __name__ == "__main__":
    main()
