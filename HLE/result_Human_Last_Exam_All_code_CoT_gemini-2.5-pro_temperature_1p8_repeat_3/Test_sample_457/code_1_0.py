import numpy as np
from scipy import stats

def solve_sba15_isotherm():
    """
    Analyzes N2 adsorption-desorption data for an SBA-15 sample to determine
    Specific Surface Area (SSA), monolayer capacity (Vm), and pore diameter (d).
    """

    # --- 1. Data Setup ---
    print("Step 1: Loading and processing the isotherm data.")

    # Provided raw data for relative pressure (P/P0) and adsorbed volume (V)
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
        463.6461, 583.7738, 466.7807, 460.9017, 455.4682, 452.4093, 449.8886,
        447.7461, 446.0029, 444.6366, 443.3456, 442.1168, 441.0445, 439.919, 438.8695,
        437.7167, 436.6012, 435.6349, 434.5094, 433.1968, 431.7222, 430.2362,
        428.62, 422.1941, 282.1162, 256.5904, 240.3565, 233.1687, 227.0982,
        224.9568, 219.4109
    ])

    # The data follows an adsorption branch (increasing pressure) then a desorption branch (decreasing pressure).
    # The split occurs at the point of maximum pressure.
    split_index = np.argmax(p_po_raw) + 1
    p_ads = p_po_raw[:split_index]
    v_ads = v_ads_raw[:split_index]
    p_des = p_po_raw[split_index:]
    v_des = v_ads_raw[split_index:]
    print("Data successfully separated into adsorption and desorption branches.\n")

    # --- 2. BET Analysis for SSA and Vm ---
    print("Step 2: Calculating Specific Surface Area (SSA) and Monolayer Capacity (Vm) using the BET method.")
    
    # Filter data for the standard BET relative pressure range (0.05 to 0.35)
    bet_mask = (p_ads >= 0.05) & (p_ads <= 0.35)
    p_bet = p_ads[bet_mask]
    v_bet = v_ads[bet_mask]

    # Calculate the BET term y = (P/P0) / [V * (1 - P/P0)]
    bet_y = p_bet / (v_bet * (1 - p_bet))
    
    # Perform linear regression on the BET plot: bet_y vs. p_bet
    slope, intercept, r_value, _, _ = stats.linregress(p_bet, bet_y)
    
    print(f"The BET equation is of the form y = mx + c.")
    print(f"Linear fit results: slope (m) = {slope:.6f}, intercept (c) = {intercept:.6f}, R^2 = {r_value**2:.5f}")

    # Calculate Vm from slope and intercept
    # Vm = 1 / (slope + intercept)
    Vm = 1 / (slope + intercept)
    print(f"Calculating monolayer capacity (Vm) using: Vm = 1 / (slope + intercept)")
    print(f"Vm = 1 / ({slope:.6f} + {intercept:.6f}) = {Vm:.1f} cm^3/g")
    
    # Calculate SSA from Vm
    # Constants for SSA calculation for N2 at 77K:
    # Avogadro's number (N_A) = 6.022e23 molecules/mol
    # Cross-sectional area of N2 (sigma) = 0.162 nm^2 = 1.62e-19 m^2
    # Molar volume of gas at STP (V_molar) = 22414 cm^3/mol
    # SSA [m^2/g] = (Vm * N_A * sigma) / V_molar = Vm * 4.353
    ssa_conversion_factor = 4.353
    SSA = Vm * ssa_conversion_factor
    print(f"Calculating Specific Surface Area (SSA) using: SSA = Vm * conversion_factor")
    print(f"SSA = {Vm:.1f} * {ssa_conversion_factor} = {SSA:.1f} m^2/g\n")

    # --- 3. Pore Diameter Analysis ---
    print("Step 3: Calculating the pore diameter (d) with the highest differential pore volume.")
    
    # Sort desorption branch by increasing pressure for analysis
    sort_indices = np.argsort(p_des)
    p_des_sorted = p_des[sort_indices]
    v_des_sorted = v_des[sort_indices]

    def calculate_diameter_nm(p_rel):
        """Calculates pore diameter in nm from relative pressure using Kelvin and Halsey equations."""
        with np.errstate(divide='ignore'):
            log_p_rel = np.log(p_rel)
            # Halsey equation for film thickness 't'
            # t[nm] = 0.354 * [5 / -ln(P/P0)]^(1/3)
            t = 0.354 * np.power(5 / -log_p_rel, 1/3)
            # Kelvin equation for pore core radius 'rk' (cylindrical pore)
            # rk[nm] = -0.9576 / ln(P/P0)
            rk = -0.9576 / log_p_rel
        return 2 * (rk + t)

    print("Using the Kelvin and Halsey equations to find pore diameter (d) from relative pressure.")
    # Calculate pore diameter for each desorption point
    d_des = calculate_diameter_nm(p_des_sorted)

    # Calculate the differential volume (dV/dD)
    dV = np.diff(v_des_sorted)
    dD = np.diff(d_des)
    
    with np.errstate(divide='ignore', invalid='ignore'):
        differential_volume = dV / dD

    # Calculate the average diameter for each interval
    d_avg = (d_des[:-1] + d_des[1:]) / 2

    # Find the diameter corresponding to the peak of the pore size distribution
    # We must ignore non-physical values (NaN, infinity) that can arise from calculations.
    valid_mask = np.isfinite(differential_volume)
    peak_index = np.argmax(differential_volume[valid_mask])
    peak_diameter = d_avg[valid_mask][peak_index]
    
    # Round to the nearest nanometer as requested
    d_rounded = int(round(peak_diameter))

    print(f"The analysis of the desorption branch shows a peak in the differential pore volume.")
    print(f"The pore diameter corresponding to this peak is {peak_diameter:.2f} nm.")
    print(f"Rounded to the nearest nanometer, the most frequent pore diameter is {d_rounded} nm.\n")
    
    # --- 4. Final Answer ---
    Vm_rounded = round(Vm, 1)
    SSA_rounded = round(SSA, 1)

    print("--- Summary of Results ---")
    print(f"Specific Surface Area (SSA): {SSA_rounded} m^2/g")
    print(f"Monolayer Capacity (Vm): {Vm_rounded} cm^3/g")
    print(f"Most Frequent Pore Diameter (d): {d_rounded} nm")
    
    final_answer_set = (SSA_rounded, Vm_rounded, d_rounded)
    print(f"\nThe final answer as an ordered set (SSA, Vm, d) is: {final_answer_set}")
    print(f"<<<{final_answer_set}>>>")


solve_sba15_isotherm()