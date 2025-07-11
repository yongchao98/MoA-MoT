import numpy as np
from scipy import stats

def solve_material_properties():
    """
    Calculates SSA, Vm, and pore diameter from N2 adsorption-desorption data.
    """
    # --- Data Preparation ---
    # Store the provided data in numpy arrays
    relative_pressure = np.array([
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

    adsorbed_volume = np.array([
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

    # Find the turning point to split adsorption and desorption branches
    turn_point_idx = np.argmax(relative_pressure)
    ads_p = relative_pressure[:turn_point_idx + 1]
    ads_v = adsorbed_volume[:turn_point_idx + 1]
    des_p = relative_pressure[turn_point_idx + 1:]
    des_v = adsorbed_volume[turn_point_idx + 1:]

    # --- Part 1 & 2: BET Specific Surface Area (SSA) and Monolayer Capacity (Vm) ---

    # Filter adsorption data for the BET relative pressure range (0.05 to 0.35)
    bet_mask = (ads_p >= 0.05) & (ads_p <= 0.35)
    bet_p = ads_p[bet_mask]
    bet_v = ads_v[bet_mask]

    # Calculate the BET plot variable y = 1 / (V * ((P0/P) - 1))
    bet_y = 1 / (bet_v * (1 / bet_p - 1))

    # Perform linear regression on the BET data
    lin_reg = stats.linregress(bet_p, bet_y)
    slope = lin_reg.slope
    intercept = lin_reg.intercept

    # Calculate monolayer capacity (Vm) from slope and intercept: Vm = 1 / (slope + intercept)
    Vm = 1 / (slope + intercept)

    # Calculate Specific Surface Area (SSA)
    # For N2 at 77K, SSA (m^2/g) = Vm (cm^3/g) * 4.353
    ssa_factor = 4.353
    SSA = Vm * ssa_factor

    # --- Part 3: Pore Diameter (d) ---

    # Reverse desorption branch to have P/P0 in increasing order for gradient calculation
    des_p_rev = des_p[::-1]
    des_v_rev = des_v[::-1]
    
    # Calculate the derivative dV/d(P/P0) to find the inflection point
    dV_dPdP0 = np.gradient(des_v_rev, des_p_rev)
    
    # Find the P/P0 corresponding to the maximum gradient (steepest part of the desorption step)
    inflection_idx = np.argmax(dV_dPdP0)
    p_inflection = des_p_rev[inflection_idx]

    # Physical constants for N2 at 77K
    gamma = 8.85e-3  # Surface tension of liquid N2 (N/m)
    V_L = 34.68e-6   # Molar volume of liquid N2 (m^3/mol)
    R = 8.314        # Ideal gas constant (J/(mol*K))
    T = 77           # Temperature (K)
    t0 = 0.354       # Thickness of a single N2 monolayer (nm)

    # Calculate Kelvin radius (r_k) using the Kelvin equation, converting m to nm
    r_k_nm = (- (2 * gamma * V_L) / (R * T * np.log(p_inflection))) * 1e9

    # Calculate the thickness of the adsorbed film (t) using the Halsey equation
    t_nm = t0 * (-5 / np.log(p_inflection))**(1/3)

    # The pore diameter is d = 2 * (r_k + t)
    d = 2 * (r_k_nm + t_nm)

    # Round d to the nearest nanometer as requested
    d_rounded = round(d)

    # --- Final Output ---
    # Report as an ordered set (SSA, Vm, d)
    print("The final answer is composed of the following values:")
    print(f"Specific Surface Area (SSA) = {SSA:.4f} m^2/g")
    print(f"Monolayer Capacity (Vm) = {Vm:.4f} cm^3/g")
    print(f"Pore Diameter (d) = {int(d_rounded)} nm")
    
    # Final answer in the specified format
    print(f"<<<({SSA:.4f}, {Vm:.4f}, {int(d_rounded)})>>>")

solve_material_properties()