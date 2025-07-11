import numpy as np

def calculate_total_voltage():
    """
    Calculates the total voltage across the load C_L considering harmonics
    and parasitic effects.
    """
    # 1. Define Constants and Parameters
    V_RF_1 = 1.0  # V, amplitude of the fundamental
    f1 = 915e6  # Hz, fundamental frequency
    Cc = 1e-12  # F
    CL = 5e-12  # F
    RL = 8e3  # Ohm
    R0_parasitic = 50.0  # Ohm
    f0_parasitic = 915e6  # Hz
    C_parasitic = 2e-15  # F

    # Effective load capacitance is the sum of C_L and C_parasitic in parallel
    C_L_eff = CL + C_parasitic

    # Harmonics to consider
    harmonics = [1, 3, 5, 7]
    V_L_rms_squared_list = []
    V_L_rms_list = []

    print("Calculating voltage contribution from each harmonic...\n")

    # 2. Loop through each harmonic
    for n in harmonics:
        # Calculate harmonic-specific values
        f_n = n * f1
        w_n = 2 * np.pi * f_n
        
        # Input voltage amplitude for the n-th harmonic
        # Voltage drops by 10% (multiplied by 0.9) for each higher harmonic step
        V_in_n = V_RF_1 * (0.9)**((n - 1) / 2)
        
        # Calculate frequency-dependent impedances (as complex numbers)
        
        # Parasitic Resistance increases with square of frequency
        R_p_n = R0_parasitic * (f_n / f0_parasitic)**2
        
        # Series coupling impedance from two C_c capacitors
        Z_c_n = 2 / (1j * w_n * Cc)
        
        # Parallel load impedance (R_L || C_L_eff)
        Z_L_n = 1 / (1/RL + 1j * w_n * C_L_eff)
        
        # Total impedance for the voltage divider calculation
        Z_total_n = R_p_n + Z_c_n + Z_L_n
        
        # Calculate voltage amplitude across the load using the voltage divider rule
        V_L_n_amp = V_in_n * abs(Z_L_n) / abs(Z_total_n)
        
        # Convert amplitude to RMS voltage
        V_L_n_rms = V_L_n_amp / np.sqrt(2)
        
        # Store results for the final calculation
        V_L_rms_list.append(V_L_n_rms)
        V_L_rms_squared_list.append(V_L_n_rms**2)
        
        print(f"--- Harmonic n={n} ---")
        print(f"  RMS Voltage at Load: V_L_{n},rms = {V_L_n_rms:.4f} V")

    # 3. Calculate the total RMS voltage
    total_V_L_rms = np.sqrt(sum(V_L_rms_squared_list))

    # 4. Print the final equation and the result
    print("\n-------------------------------------------------")
    print("The total voltage is the RMS sum of the individual harmonic voltages.")
    
    equation_str = "V_total = sqrt("
    term_str_list = [f"{v:.4f}^2" for v in V_L_rms_list]
    equation_str += " + ".join(term_str_list)
    equation_str += ")"
    
    print("\nFinal Calculation:")
    print(equation_str)
    
    sum_of_squares = sum(V_L_rms_squared_list)
    sum_of_squares_terms = " + ".join([f"{v_sq:.6f}" for v_sq in V_L_rms_squared_list])

    print(f"V_total = sqrt({sum_of_squares_terms})")
    print(f"V_total = sqrt({sum_of_squares:.6f})")
    print(f"\nFinal Total Voltage across C_L: {total_V_L_rms:.4f} V")
    print("-------------------------------------------------")

    return total_V_L_rms

if __name__ == '__main__':
    final_voltage = calculate_total_voltage()
    # The final answer is formatted below as requested.
    # print(f'<<<{final_voltage:.4f}>>>')

calculate_total_voltage()