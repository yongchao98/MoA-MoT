import numpy as np
from scipy.stats import linregress

def solve_sba15_isotherm():
    """
    Analyzes N2 adsorption-desorption data for an SBA-15 sample to determine
    specific surface area (SSA), monolayer capacity (Vm), and pore diameter (d).
    """
    # Provided data
    p_po = np.array([
        0.000097778, 0.000490286, 0.001028258, 0.005256671, 0.010218832, 0.053231858,
        0.098896253, 0.14976513, 0.200004704, 0.233517997, 0.266801197, 0.301984569,
        0.325449461, 0.35054002, 0.376514193, 0.400208966, 0.426400846, 0.450079581,
        0.476335676, 0.499897655, 0.525467416, 0.551884881, 0.575255893, 0.600463037,
        0.625966848, 0.651702807, 0.676434592, 0.701749025, 0.726670645, 0.751424402,
        0.776687236, 0.801317436, 0.826598164, 0.85144159, 0.876516006, 0.901802511,
        0.926576152, 0.96028434, 0.979849853, 0.981764776, 0.989523087, 0.999369551,
        0.984607125, 0.967573883, 0.944626529, 0.927201724, 0.90187396, 0.877154076,
        0.851603484, 0.826517066, 0.802151817, 0.777231389, 0.751936033, 0.727508526,
        0.702376073, 0.676883722, 0.650369795, 0.62669699, 0.601288503, 0.576203003,
        0.550292959, 0.525392205, 0.501925833, 0.475695455, 0.447680802, 0.423654664,
        0.382079135, 0.356158655, 0.331104613, 0.325044783, 0.299552276
    ])
    v_ads = np.array([
        34.912, 52.8679, 63.2559, 90.9391, 104.0776, 144.7389, 165.4765, 182.3897,
        196.5799, 205.2201, 213.3343, 221.7435, 227.2067, 233.1827, 239.264, 244.933,
        251.3238, 257.1772, 263.9726, 270.2935, 277.4911, 285.6177, 293.3731, 302.0746,
        312.2637, 323.6833, 335.6003, 348.5889, 362.0965, 376.0429, 390.2814, 404.5596,
        417.7411, 430.0732, 440.2639, 446.8344, 450.299, 454.6433, 458.9642, 460.8415,
        463.6461, 583.7738, 466.7807, 460.9017, 455.4682, 452.4093, 449.8886, 447.7461,
        446.0029, 444.6366, 443.3456, 442.1168, 441.0445, 439.919, 438.8695, 437.7167,
        436.6012, 435.6349, 434.5094, 433.1968, 431.7222, 430.2362, 428.62, 422.1941,
        282.1162, 256.5904, 240.3565, 233.1687, 227.0982, 224.9568, 219.4109
    ])

    # --- Part 1: BET Analysis ---
    print("--- BET Analysis ---")
    
    # Split data into adsorption (increasing pressure) and desorption branches
    # The turning point is at index 41 (P/P0 = 0.999...)
    p_ads_branch = p_po[:42]
    v_ads_branch = v_ads[:42]
    
    # Filter adsorption data for the BET relative pressure range (0.05 to 0.35)
    bet_mask = (p_ads_branch >= 0.05) & (p_ads_branch <= 0.35)
    p_bet = p_ads_branch[bet_mask]
    v_bet = v_ads_branch[bet_mask]

    # Calculate the BET transform: y = (P/P0) / (V * (1 - P/P0))
    y_bet = p_bet / (v_bet * (1 - p_bet))

    # Perform linear regression
    slope, intercept, r_value, _, _ = linregress(p_bet, y_bet)
    print(f"Linear regression on BET data (y = slope*x + intercept):")
    print(f"  y = P/[(V*(1-P/P0))], x = P/P0")
    print(f"  Slope = {slope:.6f}, Intercept = {intercept:.6f}, R^2 = {r_value**2:.4f}")

    # Calculate Vm and SSA
    # Vm = 1 / (slope + intercept)
    # SSA = Vm * N_A * A_m / V_0 = Vm * 4.353
    Vm = 1 / (slope + intercept)
    SSA = Vm * 4.353
    print("\nCalculation of Monolayer Capacity (Vm) and Specific Surface Area (SSA):")
    print(f"  Vm = 1 / (slope + intercept) = 1 / ({slope:.6f} + {intercept:.6f}) = {Vm:.2f} cm^3/g")
    print(f"  SSA = Vm * 4.353 = {Vm:.2f} * 4.353 = {SSA:.2f} m^2/g")

    # --- Part 2: Pore Diameter Analysis ---
    print("\n--- Pore Diameter Analysis (Simplified BJH) ---")

    # Constants for N2 at 77K
    T = 77.0       # Temperature (K)
    gamma = 8.85e-3  # Surface tension of liquid N2 (N/m)
    Vl = 34.68e-6  # Molar volume of liquid N2 (m^3/mol)
    R = 8.314      # Ideal gas constant (J/mol*K)

    def halsey_thickness(p_rel):
        # Halsey equation for film thickness (nm)
        # Clip p_rel to avoid log(0) or domain errors for negative logs
        p_rel = np.clip(p_rel, 1e-9, 0.999999999)
        return 0.354 * (-5 / np.log(p_rel))**(1/3.0)

    def kelvin_radius(p_rel):
        # Kelvin equation for pore core radius (nm)
        p_rel = np.clip(p_rel, 1e-9, 0.999999999)
        return -(2 * gamma * Vl * 1e9) / (R * T * np.log(p_rel))

    # Use desorption branch, sorted by increasing pressure
    p_des_branch = p_po[42:]
    v_des_branch = v_ads[42:]
    p_des_inc = p_des_branch[::-1]
    v_des_inc = v_des_branch[::-1]

    # Calculate pore diameter for each point on the desorption curve
    t_film = halsey_thickness(p_des_inc)
    r_core = kelvin_radius(p_des_inc)
    pore_diameter = 2 * (r_core + t_film)

    # Calculate differential volume dV/dD to find the peak
    dV = np.diff(v_des_inc)
    dD = np.diff(pore_diameter)
    
    # Avoid division by zero or negative dD
    valid_mask = dD > 1e-9
    dV_dD = np.zeros_like(dD)
    dV_dD[valid_mask] = dV[valid_mask] / dD[valid_mask]

    # Find the diameter corresponding to the maximum dV/dD (peak of the distribution)
    max_idx = np.argmax(dV_dD)
    # Average the diameter at the start and end of the interval with the max derivative
    d_peak = (pore_diameter[max_idx] + pore_diameter[max_idx + 1]) / 2.0
    d_rounded = int(round(d_peak))
    
    print("Calculation of Pore Diameter (d):")
    print(f"  Peak of differential pore volume (dV/dD) occurs at d = {d_peak:.2f} nm")
    print(f"  Pore diameter rounded to the nearest nanometer: {d_rounded} nm")

    # --- Final Result ---
    print("\n--- Final Result ---")
    final_answer = (round(SSA, 2), round(Vm, 2), d_rounded)
    print("Ordered set (SSA, Vm, d):")
    print(final_answer)
    return final_answer

final_tuple = solve_sba15_isotherm()
print(f'<<<{final_tuple}>>>')
