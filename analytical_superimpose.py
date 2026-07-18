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

def get_bracket(Q2, mchi):
    m_mu = 0.105658
    m_tau = 1.77686
    
    val = 1.0 # 2 neutrinos (nu_mu, nu_tau), each contributes 1/2 to the sum
    
    if Q2 > 4.0 * m_mu**2:
        val += (1.0 + 2.0 * m_mu**2 / Q2) * np.sqrt(1.0 - 4.0 * m_mu**2 / Q2)
        
    if Q2 > 4.0 * m_tau**2:
        val += (1.0 + 2.0 * m_tau**2 / Q2) * np.sqrt(1.0 - 4.0 * m_tau**2 / Q2)
        
    return val

def get_width_factor(mZp, mchi):
    val = get_bracket(mZp**2, mchi)
    
    if mZp > 2.0 * mchi:
        val += (1.0 + 2.0 * mchi**2 / mZp**2) * np.sqrt(1.0 - 4.0 * mchi**2 / mZp**2)
        
    return (mZp / (12.0 * np.pi)) * val

def sigmav_SM_factor(s, mZp, mchi, gamma_factor, gp):
    width = gp**2 * gamma_factor
    denom = (s - mZp**2)**2 + mZp**2 * width**2
    
    sm_bracket = get_bracket(s, mchi)
    initial_factor = 1.0 + 2.0 * mchi**2 / s
    
    sigmav = (s / (6.0 * np.pi * denom)) * initial_factor * sm_bracket
    return sigmav

def sigmav_ZpZp_factor(mchi, mZp):
    if mchi > mZp:
        return (1.0 / (4.0 * np.pi * mchi)) * (mchi**2 - mZp**2)**1.5 / (2.0 * mchi**2 - mZp**2)**2
    return 0.0

def get_nwa_contribution(mchi, mZp, gp, x=20.0):
    if mZp <= 2.0 * mchi:
        return 0.0
    
    gamma_factor = get_width_factor(mZp, mchi)
    Gamma_Zp = gp**2 * gamma_factor
    
    sm_bracket = get_bracket(mZp**2, mchi)
    initial_factor = 1.0 + 2.0 * mchi**2 / mZp**2
    A_val = gp**4 * (mZp**2 / (6.0 * np.pi)) * initial_factor * sm_bracket
    
    v_res_sq = 4.0 * (1.0 - 4.0 * mchi**2 / mZp**2)
    if v_res_sq <= 0:
        return 0.0
    v_res = np.sqrt(v_res_sq)
    
    bw_int = np.pi / (mZp * Gamma_Zp)
    
    dv_ds = 8.0 * mchi**2 / (mZp**4 * v_res)
    fs = (x**1.5 / (2.0 * np.sqrt(np.pi))) * v_res_sq * np.exp(-x * v_res_sq / 4.0) * dv_ds
    
    return A_val * fs * bw_int

def get_thermal_averaged_sigmav(mchi, mZp, gp, x=20.0):
    gamma_factor = get_width_factor(mZp, mchi)
    
    # We will use numerical integration for the off-resonance continuum
    # and add the exact analytical Narrow Width Approximation (NWA) for the pole
    
    v_pts = np.linspace(1e-4, 0.99, 150)
        
    integrand = []
    for v in v_pts:
        s = 4.0 * mchi**2 / (1.0 - 0.25 * v**2)
        sm_fac = sigmav_SM_factor(s, mZp, mchi, gamma_factor, gp)
        prob = (x**1.5 / (2.0 * np.sqrt(np.pi))) * v**2 * np.exp(-x * v**2 / 4.0)
        integrand.append(sm_fac * prob)
        
    sum_sm = np.trapezoid(integrand, v_pts)

    nwa = get_nwa_contribution(mchi, mZp, gp, x)
    zp_fac = sigmav_ZpZp_factor(mchi, mZp)
    
    # Note: sum_sm inherently misses the incredibly narrow peak unless v_pts happens to hit it,
    # so adding NWA is physically correct and very accurate.
    total_sigmav = gp**4 * (sum_sm + zp_fac) + nwa
    return total_sigmav

def get_relic_density(mchi, mZp, gp):
    x_f = 20.0
    T_f = mchi / x_f
    gstar = get_gstar(T_f)
    
    sigmav = get_thermal_averaged_sigmav(mchi, mZp, gp, x_f)
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
        return None

def main():
    mchi_targets = [0.05, 5.0, 500.0]
    colors = {0.05: 'red', 5.0: 'deepskyblue', 500.0: 'darkviolet'}
    
    # Reference digitized data points from Figure 2 of arXiv:2501.08622
    ref_data = {
        0.05: np.array([
            [0.001, 5e-3], [0.01, 5e-3], [0.05, 4.5e-3], [0.08, 1e-3], [0.09, 2e-4], 
            [0.095, 1e-5], [0.098, 8e-7], [0.101, 1e-5], [0.11, 2e-4], [0.12, 1e-3], 
            [0.15, 5e-3], [0.2, 8e-3], [0.5, 2.5e-2], [1.0, 6e-2], [2.0, 1.2e-1], [5.0, 3e-1]
        ]),
        5.0: np.array([
            [0.001, 3.2e-2], [0.01, 3.2e-2], [0.1, 3.2e-2], [1.0, 3.2e-2], [5.0, 2.8e-2], 
            [8.0, 8e-3], [9.0, 2e-3], [9.5, 3e-4], [9.8, 1.3e-4], [10.1, 4e-4], 
            [11.0, 3e-3], [12.0, 1e-2], [15.0, 3e-2], [30.0, 8e-2], [100.0, 3e-1]
        ]),
        500.0: np.array([
            [0.001, 3.5e-1], [0.1, 3.5e-1], [10.0, 3.5e-1], [100.0, 3.5e-1], [500.0, 3.0e-1], 
            [800.0, 8e-2], [900.0, 2.5e-2], [950.0, 8e-3], [980.0, 5.0e-3], [1010.0, 1e-2], 
            [1100.0, 6e-2], [1500.0, 2.5e-1]
        ])
    }

    plt.figure(figsize=(9, 6.5))
    ax = plt.gca()
    ax.set_facecolor('#f4f4f4')
    
    for mchi in mchi_targets:
        print(f"Scanning for mchi = {mchi} GeV...")
        mZp_min, mZp_max = 0.001, 100000.0
        
        mZp_grid = np.logspace(np.log10(mZp_min), np.log10(mZp_max), 80)
        res_center = 2.0 * mchi
        res_points = np.logspace(np.log10(res_center * 0.7), np.log10(res_center * 1.3), 40)
        mZp_grid = np.unique(np.sort(np.concatenate([mZp_grid, res_points])))
        
        results = []
        for mZp in mZp_grid:
            gp = find_gp_for_relic(mchi, mZp)
            if gp is not None:
                results.append((mZp, gp))
                
        results = np.array(results)
        
        # Plot our analytical result
        plt.plot(results[:, 0], results[:, 1], label=f"Our Calculation ($m_\\chi = {mchi}$ GeV)", color=colors[mchi], lw=2.5)
        
        # Plot paper results as dashed lines with markers
        ref = ref_data[mchi]
        plt.plot(ref[:, 0], ref[:, 1], color=colors[mchi], ls="--", marker="o", markersize=4, alpha=0.7, label=f"Paper Fig 2 ($m_\\chi = {mchi}$ GeV)")

    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(0.001, 100000.0)
    plt.ylim(1e-8, 1.0)
    plt.xlabel(r"$m_{Z'}$ (GeV)", fontsize=13)
    plt.ylabel(r"$g'$", fontsize=13)
    plt.title(r"$U(1)_{L_\mu - L_\tau}$ Relic Density: Superimposed Comparison", fontsize=13, pad=15)
    plt.grid(True, which="both", ls="--", alpha=0.6, color="white", lw=1.2)
    plt.legend(loc="lower left", fontsize=9, framealpha=0.9, ncol=2)
    
    plt.tight_layout()
    plt.savefig("relic_density_scan_superimposed.png", dpi=300)
    print("Generated relic_density_scan_superimposed.png successfully.")

if __name__ == '__main__':
    main()
