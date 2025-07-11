import numpy as np
from scipy.stats import linregress

def solve_isotherm_data():
    """
    Analyzes N2 adsorption-desorption isotherm data for an SBA-15 sample
    to calculate SSA, Vm, and pore diameter.
    """
    # 1. Data and Constants
    # Provided isotherm data
    p_po_all = np.array([
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
    v_ads_all = np.array([
        34.912, 52.8679, 63.2559, 90.9391, 104.0776, 144.7389, 165.4765, 182.3897,
        196.5799, 205.2201, 213.3343, 221.7435, 227.2067, 233.1827, 239.264,
        244.933, 251.3238, 257.1772, 263.9726, 270.2935, 277.4911, 285.6177,
        293.3731, 302.0746, 312.2637, 323.6833, 335.6003, 348.5889, 362.0965,
        376.0429, 390.2814, 404.5596, 417.7411, 430.0732, 440.2639, 446.8344,
        450.299, 454.6433, 458.9642, 460.8415, 463.6461, 583.7738, 466.7807,
        460.9017, 455.4682, 452.4093, 449.8886, 447.7461, 446.0029, 444.6366,
        443.3456, 442.1168, 441.0445, 439.919, 438.8695, 437.7167, 436.6012,
        435.6349, 434.5094, 433.1968, 431.7222, 430.2362, 428.62, 422.1941,
        282.1162, 256.5904, 240.3565, 233.1687, 227.0982, 224.9568, 219.4109
    ])

    # Physical constants for N2 at 77K
    N_A = 6.022e23         # Avogadro's number (mol^-1)
    sigma_N2 = 1.62e-19    # Cross-sectional area of N2 (m^2)
    V_molar_STP = 22414    # Molar volume of ideal gas at STP (cm^3/mol)
    GAMMA_N2 = 8.85e-3     # Surface tension of liquid N2 at 77K (N/m)
    V_L_N2 = 3.468e-5      # Molar volume of liquid N2 at 77K (m^3/mol)
    R = 8.314              # Ideal gas constant (J/(mol*K))
    T = 77                 # Temperature (K)

    # 2. Separate and sort adsorption data
    max_p_idx = np.argmax(p_po_all)
    ads_p_raw = p_po_all[:max_p_idx + 1]
    ads_v_raw = v_ads_all[:max_p_idx + 1]
    ads_sort_indices = np.argsort(ads_p_raw)
    ads_p = ads_p_raw[ads_sort_indices]
    ads_v = ads_v_raw[ads_sort_indices]

    # 3. BET Analysis
    bet_mask = (ads_p >= 0.05) & (ads_p <= 0.35)
    bet_p = ads_p[bet_mask]
    bet_v = ads_v[bet_mask]

    bet_y = bet_p / (bet_v * (1 - bet_p))
    bet_x = bet_p
    
    slope, intercept, _, _, _ = linregress(bet_x, bet_y)
    
    vm_val = 1 / (slope + intercept)
    ssa_val = (vm_val * N_A * sigma_N2) / V_molar_STP

    # 4. Pore Diameter Calculation from Adsorption Branch
    condensation_mask = (ads_p > 0.6) & (ads_p < 0.85)
    p_cond = ads_p[condensation_mask]
    v_cond = ads_v[condensation_mask]

    slopes = np.diff(v_cond) / np.diff(p_cond)
    max_slope_idx = np.argmax(slopes)
    
    p_step = (p_cond[max_slope_idx] + p_cond[max_slope_idx + 1]) / 2

    r_k_m = -(2 * GAMMA_N2 * V_L_N2) / (R * T * np.log(p_step))
    t_A = (13.99 / (0.034 - np.log10(p_step)))**0.5
    t_m = t_A * 1e-10
    
    r_p_m = r_k_m + t_m
    d_m = 2 * r_p_m
    d_nm = d_m * 1e9
    
    # 5. Final Output
    d_final = round(d_nm)
    
    # "Remember in the final code you still need to output each number in the final equation!"
    # This is interpreted as printing the values that constitute the final answer tuple.
    print(f"Monolayer Capacity (Vm) = {vm_val:.1f} cm^3/g")
    print(f"Specific Surface Area (SSA) = {ssa_val:.1f} m^2/g")
    print(f"Pore Diameter (d) = {d_nm:.1f} nm (rounded to {d_final:.0f} nm)")

    final_answer_tuple = (round(ssa_val, 1), round(vm_val, 1), d_final)
    print("\nFinal ordered set (SSA, Vm, d):")
    print(final_answer_tuple)
    
    return final_answer_tuple

# Execute the function to get the results
final_answer = solve_isotherm_data()
