import numpy as np
from scipy import stats

def analyze_isotherm_data():
    """
    Analyzes adsorption-desorption isotherm data to calculate SSA, Vm, and pore diameter.
    """
    # 1. Data Preparation
    # Provided raw data
    p_po_raw = np.array([
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
    v_ads_raw = np.array([
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

    # Find turnover point to separate adsorption and desorption branches
    turnover_idx = np.argmax(p_po_raw) + 1
    p_ads = p_po_raw[:turnover_idx]
    v_ads = v_ads_raw[:turnover_idx]
    
    # Sort desorption data by descending pressure
    p_des_raw = p_po_raw[turnover_idx:]
    v_des_raw = v_ads_raw[turnover_idx:]
    sort_indices = np.argsort(p_des_raw)[::-1]
    p_des = p_des_raw[sort_indices]
    v_des = v_des_raw[sort_indices]

    # 2. BET Analysis
    bet_mask = (p_ads >= 0.05) & (p_ads <= 0.35)
    p_bet = p_ads[bet_mask]
    v_bet = v_ads[bet_mask]

    bet_y = p_bet / (v_bet * (1 - p_bet))
    bet_x = p_bet
    slope, intercept, r_value, _, _ = stats.linregress(bet_x, bet_y)

    vm = 1 / (slope + intercept)
    c_const = (slope / intercept) + 1
    ssa = vm * 4.353  # Conversion factor for N2 at 77K

    print("--- BET Analysis ---")
    print("The linearized BET equation is: y = m*x + c")
    print("y = (P/P0) / (V * (1 - P/P0)), x = P/P0")
    print(f"Linear regression of the BET plot gives:")
    print(f"Slope (m) = (C-1)/(Vm*C) = {slope:.6f}")
    print(f"Intercept (c) = 1/(Vm*C) = {intercept:.8f}")
    print(f"Correlation coefficient (R^2) = {r_value**2:.4f}")
    print("\nSolving for Vm and C:")
    print(f"Monolayer capacity (Vm) = 1 / (slope + intercept) = 1 / ({slope:.6f} + {intercept:.8f}) = {vm:.2f} cm^3/g")
    print(f"BET constant (C) = slope/intercept + 1 = {slope:.6f}/{intercept:.8f} + 1 = {c_const:.2f}")
    print(f"\nSpecific Surface Area (SSA) = Vm * 4.353 = {vm:.2f} * 4.353 = {ssa:.2f} m^2/g")

    # 3. Pore Diameter Analysis (Simplified BJH)
    print("\n--- Pore Size Distribution Analysis (Desorption Branch) ---")
    # Constants for N2 at 77K
    T = 77; GAMMA = 8.85e-3; V_L = 34.68e-6; R = 8.314

    bjah_mask = p_des < 0.98  # Avoid instability near P/P0=1
    p_bjah = p_des[bjah_mask]
    v_bjah = v_des[bjah_mask]

    # Kelvin equation for pore core radius r_k (in nm)
    r_k = ((-2 * GAMMA * V_L) / (R * T * np.log(p_bjah))) * 1e9
    # Halsey equation for adsorbed film thickness t (in nm)
    t = 0.354 * (5 / -np.log(p_bjah))**(1/3)
    # Total pore diameter D (in nm)
    D = 2 * (r_k + t)

    # Calculate differential pore volume (dV/dD)
    dV = -np.diff(v_bjah)
    dD = -np.diff(D)
    dV_dD = dV / dD

    peak_idx = np.argmax(dV_dD)
    d_peak = (D[peak_idx] + D[peak_idx + 1]) / 2
    d_final = round(d_peak)

    print("The pore diameter is determined from the peak of the differential pore volume plot (dV/dD).")
    print(f"The maximum differential pore volume (dV/dD) is {dV_dD[peak_idx]:.2f}.")
    print(f"This peak occurs between pore diameters {D[peak_idx]:.2f} nm and {D[peak_idx+1]:.2f} nm.")
    print(f"The average pore diameter at the peak is: {d_peak:.2f} nm.")
    print(f"Rounding to the nearest nanometer gives a pore diameter (d) of: {d_final} nm.")

    # 4. Final Result
    print("\n--- Final Result ---")
    print("The final answer is an ordered set (SSA, Vm, d).")
    final_answer = (round(ssa, 2), round(vm, 2), float(d_final))
    print(f"Final Answer: {final_answer}")
    print(f"<<<{final_answer}>>>")

if __name__ == '__main__':
    analyze_isotherm_data()