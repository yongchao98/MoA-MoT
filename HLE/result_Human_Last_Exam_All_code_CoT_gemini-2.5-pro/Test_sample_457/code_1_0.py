import pandas as pd
import numpy as np
from scipy.stats import linregress

def solve_isotherm_data():
    """
    Analyzes nitrogen adsorption-desorption data to determine SSA, Vm, and pore diameter.
    """
    # --- Data Setup ---
    # The provided data is loaded into a pandas DataFrame.
    data = {
        'Relative Pressure (P/P0)': [
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
        ],
        'Adsorbed Volume (cm^3 g^-1)': [
            34.912, 52.8679, 63.2559, 90.9391, 104.0776, 144.7389, 165.4765, 182.3897,
            196.5799, 205.2201, 213.3343, 221.7435, 227.2067, 233.1827, 239.264, 244.933,
            251.3238, 257.1772, 263.9726, 270.2935, 277.4911, 285.6177, 293.3731, 302.0746,
            312.2637, 323.6833, 335.6003, 348.5889, 362.0965, 376.0429, 390.2814, 404.5596,
            417.7411, 430.0732, 440.2639, 446.8344, 450.299, 454.6433, 458.9642, 460.8415,
            463.6461, 583.7738, 466.7807, 460.9017, 455.4682, 452.4093, 449.8886, 447.7461,
            446.0029, 444.6366, 443.3456, 442.1168, 441.0445, 439.919, 438.8695, 437.7167,
            436.6012, 435.6349, 434.5094, 433.1968, 431.7222, 430.2362, 428.62, 422.1941,
            282.1162, 256.5904, 240.3565, 233.1687, 227.0982, 224.9568, 219.4109
        ]
    }
    df = pd.DataFrame(data)
    df.columns = ['p_p0', 'volume']

    # --- Step 1: Split data into adsorption and desorption branches ---
    print("Step 1: Splitting data into adsorption and desorption branches.")
    adsorption_end_index = df['p_p0'].idxmax()
    adsorption_df = df.iloc[0:adsorption_end_index + 1]
    desorption_df = df.iloc[adsorption_end_index + 1:].sort_values(by='p_p0').reset_index(drop=True)
    print(f"Adsorption branch contains {len(adsorption_df)} points.")
    print(f"Desorption branch contains {len(desorption_df)} points.\n")

    # --- Step 2: BET Analysis ---
    print("Step 2: Performing BET analysis to find Vm and SSA.")
    
    # Constants for BET
    N_A = 6.022e23  # Avogadro's number, molecules/mol
    sigma_n2 = 1.62e-19  # Cross-sectional area of N2, m^2/molecule
    V_molar_STP = 22414  # Molar volume of gas at STP, cm^3/mol
    
    # Select data in the BET range (0.05 to 0.35)
    bet_df = adsorption_df[(adsorption_df['p_p0'] >= 0.05) & (adsorption_df['p_p0'] <= 0.35)].copy()
    
    # Calculate BET variable y = 1 / (V * ((P0/P) - 1))
    bet_df['bet_y'] = bet_df['p_p0'] / (bet_df['volume'] * (1 - bet_df['p_p0']))
    
    # Perform linear regression
    slope, intercept, r_value, _, _ = linregress(bet_df['p_p0'], bet_df['bet_y'])
    
    print("Linear regression performed on BET data in the P/P0 range [0.05, 0.35].")
    
    # Calculate Vm and C constant
    Vm = 1 / (slope + intercept)
    C = slope / intercept + 1
    
    # Calculate SSA
    SSA = (Vm * N_A * sigma_n2) / V_molar_STP
    
    print("\n--- BET Results ---")
    print(f"Monolayer capacity (Vm) = {Vm:.4f} cm^3/g")
    print(f"Specific Surface Area (SSA) = {SSA:.4f} m^2/g")
    print(f"The BET equation for this data is:")
    print(f"1 / (V * ((P0/P) - 1)) = {slope:.4f} * (P/P0) + {intercept:.6f}\n")

    # --- Step 3: Pore Size Distribution (BJH-like method) ---
    print("Step 3: Performing pore size analysis on the desorption branch.")
    
    # Constants for Kelvin equation (Nitrogen at 77 K)
    kelvin_const_rk = 0.9577  # nm

    psd_df = desorption_df[desorption_df['p_p0'] < 0.98].copy()

    with np.errstate(divide='ignore'):
        psd_df['ln_p_p0'] = np.log(psd_df['p_p0'])
        psd_df['rk'] = -kelvin_const_rk / psd_df['ln_p_p0']
        psd_df['t'] = 0.354 * (-5.0 / psd_df['ln_p_p0'])**(1/3)

    psd_df['Dp'] = 2 * (psd_df['rk'] + psd_df['t'])
    
    delta_V = np.diff(psd_df['volume'])
    delta_Dp = np.diff(psd_df['Dp'])
    
    dV_dDp = np.divide(delta_V, delta_Dp, out=np.zeros_like(delta_V), where=delta_Dp!=0)
    
    avg_Dp = (psd_df['Dp'].iloc[:-1].values + psd_df['Dp'].iloc[1:].values) / 2
    
    peak_index = np.argmax(dV_dDp)
    peak_pore_diameter = avg_Dp[peak_index]

    print("Pore size distribution calculated using a simplified BJH method.")
    print(f"The pore diameter with the highest differential pore volume is {peak_pore_diameter:.4f} nm.")
    
    peak_pore_diameter_rounded = round(peak_pore_diameter)
    print(f"Rounded to the nearest nanometer, the most probable pore diameter is {peak_pore_diameter_rounded} nm.\n")

    # --- Final Answer ---
    final_ssa = round(SSA, 2)
    final_vm = round(Vm, 2)
    final_d = int(peak_pore_diameter_rounded)
    
    final_answer = (final_ssa, final_vm, final_d)
    
    print("--- Summary of Results ---")
    print(f"Final ordered set (SSA, Vm, d): {final_answer}")
    
    return final_answer

# Execute the function and capture the final answer
final_tuple_result = solve_isotherm_data()
print(f"<<<{final_tuple_result}>>>")