import numpy as np
import pandas as pd
from scipy.stats import linregress

def solve_isotherm_analysis():
    """
    Analyzes adsorption-desorption isotherm data to determine SSA, Vm, and pore diameter.
    """
    # 1. Data Preparation
    data = {
        'P/P0': [0.000097778, 0.000490286, 0.001028258, 0.005256671, 0.010218832, 0.053231858, 0.098896253, 0.14976513, 0.200004704, 0.233517997, 0.266801197, 0.301984569, 0.325449461, 0.35054002, 0.376514193, 0.400208966, 0.426400846, 0.450079581, 0.476335676, 0.499897655, 0.525467416, 0.551884881, 0.575255893, 0.600463037, 0.625966848, 0.651702807, 0.676434592, 0.701749025, 0.726670645, 0.751424402, 0.776687236, 0.801317436, 0.826598164, 0.85144159, 0.876516006, 0.901802511, 0.926576152, 0.96028434, 0.979849853, 0.981764776, 0.989523087, 0.999369551, 0.984607125, 0.967573883, 0.944626529, 0.927201724, 0.90187396, 0.877154076, 0.851603484, 0.826517066, 0.802151817, 0.777231389, 0.751936033, 0.727508526, 0.702376073, 0.676883722, 0.650369795, 0.62669699, 0.601288503, 0.576203003, 0.550292959, 0.525392205, 0.501925833, 0.475695455, 0.447680802, 0.423654664, 0.382079135, 0.356158655, 0.331104613, 0.325044783, 0.299552276],
        'V': [34.912, 52.8679, 63.2559, 90.9391, 104.0776, 144.7389, 165.4765, 182.3897, 196.5799, 205.2201, 213.3343, 221.7435, 227.2067, 233.1827, 239.264, 244.933, 251.3238, 257.1772, 263.9726, 270.2935, 277.4911, 285.6177, 293.3731, 302.0746, 312.2637, 323.6833, 335.6003, 348.5889, 362.0965, 376.0429, 390.2814, 404.5596, 417.7411, 430.0732, 440.2639, 446.8344, 450.299, 454.6433, 458.9642, 460.8415, 463.6461, 583.7738, 466.7807, 460.9017, 455.4682, 452.4093, 449.8886, 447.7461, 446.0029, 444.6366, 443.3456, 442.1168, 441.0445, 439.919, 438.8695, 437.7167, 436.6012, 435.6349, 434.5094, 433.1968, 431.7222, 430.2362, 428.62, 422.1941, 282.1162, 256.5904, 240.3565, 233.1687, 227.0982, 224.9568, 219.4109]
    }
    df = pd.DataFrame(data)

    # Find the index of the highest pressure point to split the data
    split_idx = df['P/P0'].idxmax()
    ads_df = df.iloc[:split_idx + 1]
    des_df = df.iloc[split_idx + 1:]

    # --- 2. BET Analysis for SSA and Vm ---
    print("--- BET Analysis ---")
    N_A = 6.022e23  # Avogadro's number, molecules/mol
    A_m_N2 = 0.162e-18  # Cross-sectional area of N2, m^2/molecule
    V_molar = 22414  # Molar volume of ideal gas at STP, cm^3/mol

    bet_df = ads_df[(ads_df['P/P0'] > 0.05) & (ads_df['P/P0'] < 0.35)].copy()
    bet_df['x'] = bet_df['P/P0']
    bet_df['y'] = 1 / (bet_df['V'] * (1 / bet_df['P/P0'] - 1))

    slope, intercept, r_value, _, _ = linregress(bet_df['x'], bet_df['y'])

    Vm = 1 / (slope + intercept)
    C = slope / intercept + 1
    SSA = (Vm * N_A * A_m_N2) / V_molar

    print(f"The linearized BET equation is: 1/(V*(P0/P-1)) = {slope:.5f} * (P/P0) + {intercept:.5f}")
    print(f"Correlation coefficient R^2: {r_value**2:.5f}\n")
    print(f"Monolayer Capacity (Vm) = 1 / ({slope:.5f} + {intercept:.5f}) = {Vm:.4f} cm^3/g")
    print(f"Specific Surface Area (SSA) = ({Vm:.4f} * {N_A} * {A_m_N2}) / {V_molar} = {SSA:.4f} m^2/g\n")

    # --- 3. Pore Diameter Analysis (Simplified BJH on desorption branch) ---
    print("--- Pore Diameter Analysis ---")
    des_df = des_df.sort_values(by='P/P0', ascending=False).reset_index(drop=True)
    
    # Constants for N2 at 77K
    gamma = 8.85e-3      # Surface tension (N/m)
    V_L = 3.468e-5       # Liquid molar volume (m^3/mol)
    R = 8.314            # Ideal gas constant (J/mol·K)
    T = 77               # Temperature (K)
    LDCF_cm = 34.68 / 22414 # Gas to liquid volume conversion factor

    # Calculate pore diameter for each point
    p_vals = des_df['P/P0'].values
    ln_p = np.log(p_vals)
    
    # Kelvin radius in meters
    r_k = -(2 * gamma * V_L) / (R * T * ln_p)
    # Halsey thickness of adsorbed layer in meters
    t = (0.354e-9) * ((-5 / ln_p)**(1/3))
    
    # Pore diameter in nanometers
    D_p_nm = 2 * (r_k + t) * 1e9
    des_df['D_p_nm'] = D_p_nm
    
    # Calculate differential pore volume (dV_liq / dD_p)
    V_ads_liq = des_df['V'].values * LDCF_cm
    
    dV_liq = -np.diff(V_ads_liq) # Change in desorbed liquid volume
    dD_p = -np.diff(D_p_nm)     # Change in pore diameter

    # Avoid division by zero
    dD_p[dD_p == 0] = 1e-9
    
    diff_pore_vol = dV_liq / dD_p
    
    # Associate the differential volume with the average diameter of the interval
    avg_D_p = (D_p_nm[:-1] + D_p_nm[1:]) / 2

    # Filter for the meaningful part of the curve (often < 15 nm for SBA-15)
    valid_range = (avg_D_p > 2) & (avg_D_p < 15)
    
    if np.any(valid_range):
        peak_idx = np.argmax(diff_pore_vol[valid_range])
        peak_diameter = avg_D_p[valid_range][peak_idx]
    else: # Fallback if filtering is too strict
        peak_idx = np.argmax(diff_pore_vol)
        peak_diameter = avg_D_p[peak_idx]
        
    d_rounded = round(peak_diameter)
    
    print(f"The pore diameter with the highest differential pore volume is {peak_diameter:.4f} nm.")
    print(f"Rounded to the nearest nanometer, d = {d_rounded} nm.\n")

    # --- 4. Final Output ---
    result_tuple = (SSA, Vm, float(d_rounded))
    print("--- Final Result ---")
    print("The ordered set (SSA, Vm, d) is:")
    print(result_tuple)
    
    # Final answer in the specified format for the platform
    print(f"<<<{result_tuple}>>>")

# Run the analysis
solve_isotherm_analysis()