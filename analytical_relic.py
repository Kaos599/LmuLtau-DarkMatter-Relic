import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar

def get_gstar(T):
    if T < 0.001:
        return 3.36
    elif T < 0.1:
        return 10.75
    elif T < 1.0:
        return 61.75
    else:
        return 86.25

def get_width_factor(mZp, mchi):
    m_mu = 0.105658
    m_tau = 1.77686
    
    gamma = mZp / (12.0 * np.pi)
    
    if mZp > 2.0 * m_mu:
        gamma += (mZp / (12.0 * np.pi)) * (1.0 + 2.0 * m_mu**2 / mZp**2) * np.sqrt(1.0 - 4.0 * m_mu**2 / mZp**2)
        
    if mZp > 2.0 * m_tau:
        gamma += (mZp / (12.0 * np.pi)) * (1.0 + 2.0 * m_tau**2 / mZp**2) * np.sqrt(1.0 - 4.0 * m_tau**2 / mZp**2)
        
    if mZp > 2.0 * mchi:
        gamma += (mZp / (12.0 * np.pi)) * (1.0 + 2.0 * mchi**2 / mZp**2) * np.sqrt(1.0 - 4.0 * mchi**2 / mZp**2)
        
    return gamma

def sigmav_SM_factor(s, mZp, mchi, gamma_factor, gp):
    m_mu = 0.105658
    m_tau = 1.77686
    
    width = gp**2 * gamma_factor
    denom = (s - mZp**2)**2 + mZp**2 * width**2
    
    sigmav = 0.0
    sigmav += s / (2.0 * np.pi)
    
    if s > 4.0 * m_mu**2:
        sigmav += (1.0 / (2.0 * np.pi * np.sqrt(s))) * (s + 2.0 * m_mu**2) * np.sqrt(s - 4.0 * m_mu**2)
        
    if s > 4.0 * m_tau**2:
        sigmav += (1.0 / (2.0 * np.pi * np.sqrt(s))) * (s + 2.0 * m_tau**2) * np.sqrt(s - 4.0 * m_tau**2)
        
    return sigmav / denom

def sigmav_ZpZp_factor(mchi, mZp):
    if mchi > mZp:
        return (1.0 / (4.0 * np.pi * mchi)) * (mchi**2 - mZp**2)**1.5 / (2.0 * mchi**2 - mZp**2)**2
    return 0.0

def get_thermal_averaged_sigmav(mchi, mZp, gp, x=20.0):
    gamma_factor = get_width_factor(mZp, mchi)
    
    v_res = None
    if mZp > 2.0 * mchi:
        val = 1.0 - 4.0 * mchi**2 / mZp**2
        if val >= 0:
            v_res = 2.0 * np.sqrt(val)
            
    if v_res is not None and v_res < 1.0:
        v_pts = np.concatenate([
            np.linspace(1e-4, max(1e-4, v_res - 0.04), 40),
            np.linspace(max(1e-4, v_res - 0.04), min(0.99, v_res + 0.04), 80),
            np.linspace(min(0.99, v_res + 0.04), 0.99, 40)
        ])
        v_pts = np.unique(np.sort(v_pts))
    else:
        v_pts = np.linspace(1e-4, 0.99, 100)
        
    integrand = []
    for v in v_pts:
        s = 4.0 * mchi**2 / (1.0 - 0.25 * v**2)
        sm_fac = sigmav_SM_factor(s, mZp, mchi, gamma_factor, gp)
        prob = (x**1.5 / (8.0 * np.sqrt(np.pi))) * v**2 * np.exp(-x * v**2 / 4.0)
        integrand.append(sm_fac * prob)
        
    sum_sm = np.trapezoid(integrand, v_pts)

    zp_fac = sigmav_ZpZp_factor(mchi, mZp)
    
    total_sigmav = gp**4 * (sum_sm + zp_fac)
    return total_sigmav

def get_relic_density(mchi, mZp, gp):
    x_f = 20.0
    T_f = mchi / x_f
    gstar = get_gstar(T_f)
    
    sigmav = get_thermal_averaged_sigmav(mchi, mZp, gp, x_f)
    
    # Omega h^2 = 2 * (1.07e9 / (M_Pl * sqrt(g*) * sigmav / x_f))
    # sigmav must be in GeV^-2 here to match M_Pl in GeV
    M_Pl = 1.22e19
    omega = 2.14e9 * x_f / (M_Pl * np.sqrt(gstar) * sigmav)
    return omega


def objective_function(log_gp, mchi, mZp):
    gp = 10**log_gp
    omega = get_relic_density(mchi, mZp, gp)
    return np.log10(omega) - np.log10(0.12)

def find_gp_for_relic(mchi, mZp):
    try:
        f_min = objective_function(-8.0, mchi, mZp)
        f_max = objective_function(0.3, mchi, mZp)
        
        if f_min * f_max > 0:
            if f_max > 0:
                return None
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
        print(f"Error in find_gp_for_relic for mchi={mchi}, mZp={mZp}: {e}")
        import traceback
        traceback.print_exc()
        return None


def main():
    mchi_targets = [0.05, 5.0, 500.0]
    colors = {0.05: 'red', 5.0: 'deepskyblue', 500.0: 'darkviolet'}
    
    plt.figure(figsize=(9, 6.5))
    
    # Set background style similar to the figure
    ax = plt.gca()
    ax.set_facecolor('#f4f4f4')
    
    for mchi in mchi_targets:
        print(f"Scanning for mchi = {mchi} GeV...")
        mZp_min, mZp_max = 0.001, 100000.0
        
        # Base grid
        mZp_grid = np.logspace(np.log10(mZp_min), np.log10(mZp_max), 80)
        
        # Insert resonance points dense
        res_center = 2.0 * mchi
        res_points = np.logspace(np.log10(res_center * 0.7), np.log10(res_center * 1.3), 40)
        mZp_grid = np.unique(np.sort(np.concatenate([mZp_grid, res_points])))
        
        results = []
        for mZp in mZp_grid:
            gp = find_gp_for_relic(mchi, mZp)
            if gp is not None:
                results.append((mZp, gp))
                
        results = np.array(results)
        plt.plot(results[:, 0], results[:, 1], label=f"$m_\\chi = {mchi}$ GeV", color=colors[mchi], lw=2.5)
        
        # Mark the resonance minimum with text label matching Figure 2 style
        min_idx = np.argmin(results[:, 1])
        res_mZp = results[min_idx, 0]
        res_gp = results[min_idx, 1]
        plt.text(res_mZp * 1.1, res_gp * 1.5, f"$m_\\chi={mchi}\\mathrm{{GeV}}$", color=colors[mchi], fontsize=10, weight='bold')

    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(0.001, 100000.0)
    plt.ylim(1e-8, 1.0)
    plt.xlabel(r"$m_{Z'}$ (GeV)", fontsize=13)
    plt.ylabel(r"$g'$", fontsize=13)
    plt.title(r"$U(1)_{L_\mu - L_\tau}$ Relic Density Parameter Space ($\Omega_\chi h^2 = 0.12$)", fontsize=13, pad=15)
    plt.grid(True, which="both", ls="--", alpha=0.6, color="white", lw=1.2)
    plt.legend(loc="lower right", fontsize=11, framealpha=0.9)

    
    plt.tight_layout()
    plt.savefig("relic_density_scan_analytical.png", dpi=300)
    print("Generated relic_density_scan_analytical.png successfully.")

if __name__ == '__main__':
    main()
