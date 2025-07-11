import numpy as np
from scipy import stats

def solve_isotherm_data():
    """
    Analyzes adsorption-desorption isotherm data to find SSA, Vm, and pore diameter.
    """
    # --- Data Setup ---
    P_Po_all = np.array([
        0.000097778, 0.000490286, 0.001028258, 0.005256671, 0.010218832, 0.053231858,
        0.098896253, 0.14976513,  0.200004704, 0.233517997, 0.266801197, 0.301984569,
        0.325449461, 0.35054002,  0.376514193, 0.400208966, 0.426400846, 0.450079581,
        0.476335676, 0.499897655, 0.525467416, 0.551884881, 0.575255893, 0.600463037,
        0.625966848, 0.651702807, 0.676434592, 0.701749025, 0.726670645, 0.751424402,
        0.776687236, 0.801317436, 0.826598164, 0.85144159,  0.876516006, 0.901802511,
        0.926576152, 0.96028434,  0.979849853, 0.981764776, 0.989523087, 0.999369551,
        0.984607125, 0.967573883, 0.944626529, 0.927201724, 0.90187396,  0.877154076,
        0.851603484, 0.826517066, 0.802151817, 0.777231389, 0.751936033, 0.727508526,
        0.702376073, 0.676883722, 0.650369795, 0.62669699,  0.601288503, 0.576203003,
        0.550292959, 0.525392205, 0.501925833, 0.475695455, 0.447680802, 0.423654664,
        0.382079135, 0.356158655, 0.331104613, 0.325044783, 0.299552276
    ])
    V_ads_all = np.array([
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

    # --- Part 1: BET Analysis (SSA and Vm) ---
    print("--- Part 1: BET Analysis for Specific Surface Area (SSA) and Monolayer Capacity (Vm) ---")
    max_p_idx = np.argmax(P_Po_all)
    P_ads = P_Po_all[:max_p_idx + 1]
    V_ads = V_ads_all[:max_p_idx + 1]

    bet_mask = (P_ads > 0.05) & (P_ads < 0.35)
    P_bet = P_ads[bet_mask]
    V_bet = V_ads[bet_mask]

    bet_y = P_bet / (V_bet * (1 - P_bet))
    slope, intercept, r_value, _, _ = stats.linregress(P_bet, bet_y)
    
    print("\nStep 1.1: Linear Regression using the BET equation:")
    print(f"y = m*x + c, where y = (P/P0)/[V(1-P/P0)] and x = P/P0")
    print(f"The linear fit gives: slope(m) = {slope:.5f}, intercept(c) = {intercept:.5f}")
    print(f"The goodness of fit, R-squared = {r_value**2:.4f}")

    print("\nStep 1.2: Calculate Monolayer Capacity (Vm)")
    Vm = 1 / (slope + intercept)
    print(f"Vm = 1 / (m + c) = 1 / ({slope:.5f} + {intercept:.5f}) = {Vm:.2f} cm^3/g")
    
    print("\nStep 1.3: Calculate Specific Surface Area (SSA)")
    N_A = 6.022e23  # Avogadro's number, mol^-1
    sigma_N2 = 1.62e-19  # Cross-sectional area of N2, m^2
    V0 = 22414  # Molar volume of gas at STP, cm^3/mol
    SSA = Vm * (N_A * sigma_N2 / V0)
    conversion_factor = (N_A * sigma_N2 / V0)
    print(f"SSA = Vm * (N_A * sigma_N2 / V0) = {Vm:.2f} * ({N_A:.3e} * {sigma_N2:.2e} / {V0})")
    print(f"SSA = {Vm:.2f} * {conversion_factor:.3f} = {SSA:.2f} m^2/g")
    
    # --- Part 2: Pore Diameter Analysis ---
    print("\n--- Part 2: Pore Diameter (d) Analysis ---")
    P_des = P_Po_all[max_p_idx + 1:]
    V_des = V_des_all[max_p_idx + 1:]

    dV_dP = np.diff(V_des) / np.diff(P_des)
    max_dvdp_idx = np.argmax(dV_dP)
    P_Po_peak = (P_des[max_dvdp_idx] + P_des[max_dvdp_idx+1]) / 2

    print("\nStep 2.1: Find the pressure of maximum pore emptying from the desorption branch")
    print(f"The steepest desorption occurs at P/P0 = {P_Po_peak:.4f}")

    print("\nStep 2.2: Calculate Kelvin Radius (r_k) and Adsorbed Layer Thickness (t)")
    gamma_N2 = 8.85e-3; V_L_N2 = 3.468e-5; R = 8.314; T = 77
    r_k = -(2 * gamma_N2 * V_L_N2) / (R * T * np.log(P_Po_peak))
    t_nm = 0.354 * (-5 / np.log(P_Po_peak))**(1/3)
    r_k_nm = r_k * 1e9

    print(f"Kelvin equation: r_k = -(2 * gamma * V_L) / (R * T * ln(P/P0))")
    print(f"r_k = -(2 * {gamma_N2} * {V_L_N2}) / ({R} * {T} * ln({P_Po_peak:.4f})) = {r_k_nm:.2f} nm")
    print(f"Thickness equation: t = 0.354 * [-5 / ln(P/P0)]^(1/3)")
    print(f"t = 0.354 * [-5 / ln({P_Po_peak:.4f})]^(1/3) = {t_nm:.2f} nm")

    print("\nStep 2.3: Calculate final Pore Diameter (d)")
    r_p_nm = r_k_nm + t_nm
    d_nm = 2 * r_p_nm
    d_final = round(d_nm)
    print(f"Pore Radius (r_p) = r_k + t = {r_k_nm:.2f} + {t_nm:.2f} = {r_p_nm:.2f} nm")
    print(f"Pore Diameter (d) = 2 * r_p = 2 * {r_p_nm:.2f} = {d_nm:.2f} nm")
    print(f"Rounding to the nearest nanometer, d = {d_final:.0f} nm")

    # --- Final Result ---
    print("\n--- Final Ordered Set (SSA, Vm, d) ---")
    final_tuple = (SSA, Vm, float(d_final))
    print(f"({final_tuple[0]:.2f}, {final_tuple[1]:.2f}, {final_tuple[2]:.0f})")

    # Required output format
    print(f"<<<({final_tuple[0]:.2f}, {final_tuple[1]:.2f}, {final_tuple[2]:.0f})>>>")

solve_isotherm_data()