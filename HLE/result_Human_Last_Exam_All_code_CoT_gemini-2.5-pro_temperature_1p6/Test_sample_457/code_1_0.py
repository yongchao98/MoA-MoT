import numpy as np
from scipy import stats

def solve_isotherm_data():
    """
    Analyzes nitrogen adsorption-desorption isotherm data to determine
    Specific Surface Area (SSA), Monolayer Capacity (Vm), and dominant Pore Diameter (d).
    """
    # --- Step 1: Set up the data and constants ---

    # Provided data
    p_po_full = np.array([
        0.000097778, 0.000490286, 0.001028258, 0.005256671, 0.010218832,
        0.053231858, 0.098896253, 0.14976513,  0.200004704, 0.233517997,
        0.266801197, 0.301984569, 0.325449461, 0.35054002,  0.376514193,
        0.400208966, 0.426400846, 0.450079581, 0.476335676, 0.499897655,
        0.525467416, 0.551884881, 0.575255893, 0.600463037, 0.625966848,
        0.651702807, 0.676434592, 0.701749025, 0.726670645, 0.751424402,
        0.776687236, 0.801317436, 0.826598164, 0.85144159,  0.876516006,
        0.901802511, 0.926576152, 0.96028434,  0.979849853, 0.981764776,
        0.989523087, 0.999369551, 0.984607125, 0.967573883, 0.944626529,
        0.927201724, 0.90187396,  0.877154076, 0.851603484, 0.826517066,
        0.802151817, 0.777231389, 0.751936033, 0.727508526, 0.702376073,
        0.676883722, 0.650369795, 0.62669699,  0.601288503, 0.576203003,
        0.550292959, 0.525392205, 0.501925833, 0.475695455, 0.447680802,
        0.423654664, 0.382079135, 0.356158655, 0.331104613, 0.325044783,
        0.299552276
    ])
    v_full = np.array([
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

    # Physical constants for calculations
    SSA_FACTOR = 4.353  # Conversion factor for N2 at 77K: (N_A * sigma_N2) / V_molar_gas
    GAMMA = 8.85e-3      # Surface tension of liquid N2 at 77K [N/m]
    VL = 34.68e-6      # Molar volume of liquid N2 at 77K [m^3/mol]
    R = 8.31446        # Ideal gas constant [J/(mol·K)]
    T = 77             # Temperature [K]

    # Separate data into adsorption and desorption branches
    split_idx = np.argmax(p_po_full)
    p_ads = p_po_full[:split_idx + 1]
    v_ads = v_full[:split_idx + 1]
    p_des = p_po_full[split_idx + 1:]
    v_des = v_full[split_idx + 1:]

    # --- Step 2: BET Analysis for SSA and Vm ---
    print("--- BET Analysis ---")
    bet_mask = (p_ads > 0.05) & (p_ads < 0.35)
    p_bet = p_ads[bet_mask]
    v_bet = v_ads[bet_mask]

    bet_y = p_bet / (v_bet * (1 - p_bet))
    slope, intercept, r_value, _, _ = stats.linregress(p_bet, bet_y)
    
    Vm = 1 / (slope + intercept)
    SSA = Vm * SSA_FACTOR
    
    print(f"1. Linear fit of BET equation y = slope * x + intercept, where y = (P/P0)/(V(1-P/P0)) and x = P/P0:")
    print(f"   y = {slope:.5f} * x + {intercept:.5f} (R^2 = {r_value**2:.4f})")
    print("\n2. Monolayer Capacity (Vm) Calculation:")
    print(f"   Vm = 1 / (slope + intercept) = 1 / ({slope:.5f} + {intercept:.5f}) = {Vm:.2f} cm^3/g")
    print("\n3. Specific Surface Area (SSA) Calculation:")
    print(f"   SSA = Vm * conversion_factor = {Vm:.2f} * {SSA_FACTOR} = {SSA:.2f} m^2/g")

    # --- Step 3: Pore Diameter Analysis ---
    print("\n--- Pore Diameter Analysis ---")
    dv_dp = np.diff(v_des) / np.diff(p_des)
    p_avg = (p_des[:-1] + p_des[1:]) / 2
    peak_idx = np.argmax(dv_dp)
    p_peak = p_avg[peak_idx]

    rk_numerator = -2 * GAMMA * VL * 1e9  # Factor 1e9 converts m to nm
    rk_denominator = R * T * np.log(p_peak)
    rk = rk_numerator / rk_denominator

    t = 0.354 * (5 / (-np.log(p_peak)))**(1/3)

    d = 2 * (rk + t)
    d_rounded = round(d)

    print(f"1. Pressure at peak desorption (from dV/d(P/P0) analysis): P/P0 = {p_peak:.4f}")
    print("\n2. Kelvin Radius (rk) Calculation:")
    print(f"   rk = -(2 * {GAMMA:.3e} * {VL:.3e} * 1e9) / ({R:.3f} * {T} * ln({p_peak:.4f})) = {rk:.2f} nm")
    print("\n3. Adsorbed Layer Thickness (t) Calculation (Halsey):")
    print(f"   t = 0.354 * [5 / -ln({p_peak:.4f})]^(1/3) = {t:.2f} nm")
    print("\n4. Pore Diameter (d) Calculation:")
    print(f"   d = 2 * (rk + t) = 2 * ({rk:.2f} + {t:.2f}) = {d:.2f} nm")
    print(f"   Rounded pore diameter = {d_rounded} nm")
    
    # --- Step 4: Final Answer ---
    final_answer = (round(SSA, 2), round(Vm, 2), int(d_rounded))
    print("\n--- Final Ordered Set (SSA, Vm, d) ---")
    print(final_answer)
    print(f"<<<{final_answer}>>>")

solve_isotherm_data()