import numpy as np

def solve_isotherm_data():
    """
    Analyzes N2 adsorption-desorption isotherm data for an SBA-15 sample
    to calculate SSA, Vm, and pore diameter.
    """
    # Provided data
    p_po_full = np.array([
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

    V_full = np.array([
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
    
    # Split data into adsorption and desorption branches
    # Adsorption is up to the highest P/P0 value (index 41)
    p_ads = p_po_full[:42]
    V_ads = V_full[:42]
    # Desorption is the rest of the data
    p_des = p_po_full[42:]
    V_des = V_full[42:]

    # --- Part 1: BET Analysis ---
    print("--- BET Surface Area Analysis ---")
    
    # Filter for BET range (0.05 <= P/P0 <= 0.35)
    bet_mask = (p_ads >= 0.05) & (p_ads <= 0.35)
    p_bet = p_ads[bet_mask]
    V_bet = V_ads[bet_mask]

    # BET transformation: y = 1 / (V * (P0/P - 1)) , x = P/P0
    bet_y = 1 / (V_bet * (1/p_bet - 1))
    
    # Linear regression
    slope, intercept = np.polyfit(p_bet, bet_y, 1)

    # Calculate Vm and C constant
    Vm = 1 / (slope + intercept)
    C = slope / intercept + 1

    # Calculate SSA
    N_A = 6.022e23  # Avogadro's number (mol^-1)
    a_m = 0.162e-18 # Cross-sectional area of N2 (m^2)
    V_molar = 22414 # Molar volume of gas at STP (cm^3/mol)
    SSA = (Vm * N_A * a_m) / V_molar

    print("The linear BET equation is: 1/[V(P0/P-1)] = {:.6f} * (P/P0) + {:.6f}".format(slope, intercept))
    print("Monolayer Capacity (Vm): {:.2f} cm^3/g".format(Vm))
    print("Specific Surface Area (SSA): {:.2f} m^2/g".format(SSA))
    print("BET Constant (C): {:.2f}".format(C))
    
    # --- Part 2: Pore Diameter Analysis ---
    print("\n--- Pore Size Analysis (BJH Method on Desorption Branch) ---")
    
    # Sort desorption data by P/P0 for derivative calculation
    sort_indices = np.argsort(p_des)
    p_des_sorted = p_des[sort_indices]
    V_des_sorted = V_des[sort_indices]
    
    # Calculate derivative dV/d(P/P0)
    dV_dp = np.gradient(V_des_sorted, p_des_sorted)
    
    # Find the P/P0 corresponding to the peak of the derivative (inflection point)
    peak_idx = np.argmax(dV_dp)
    p_peak = p_des_sorted[peak_idx]

    # Calculate pore diameter using BJH equations for cylindrical pores
    ln_p = np.log(p_peak)
    r_k = -0.953 / ln_p  # Kelvin radius in nm
    t = 0.354 * (-5 / ln_p)**(1/3) # Adsorbed layer thickness in nm
    d_p = 2 * (r_k + t) # Pore diameter in nm
    d_p_rounded = int(round(d_p))
    
    print("The peak of the differential pore volume occurs at P/P0 = {:.4f}".format(p_peak))
    print("Calculated Pore Diameter (d): {:.2f} nm".format(d_p))
    print("Pore Diameter rounded to the nearest nanometer: {} nm".format(d_p_rounded))
    
    # --- Final Answer ---
    # Report as an ordered set (SSA, Vm, d)
    final_answer = (round(SSA, 2), round(Vm, 2), d_p_rounded)
    print("\nFinal ordered set (SSA, Vm, d):")
    # Using a f-string to make the output clear.
    print(f"({final_answer[0]}, {final_answer[1]}, {final_answer[2]})")

    return final_answer

if __name__ == '__main__':
    result = solve_isotherm_data()
    # The final answer format for the platform
    print(f"\n<<<{result}>>>")