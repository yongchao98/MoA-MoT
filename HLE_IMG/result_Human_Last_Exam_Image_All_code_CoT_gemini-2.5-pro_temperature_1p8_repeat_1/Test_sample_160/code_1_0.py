import numpy as np

def calculate_total_voltage():
    """
    Calculates the total voltage across the load capacitor CL, considering
    harmonic distortions and parasitic effects.
    """
    # --- 1. Parse Parameters ---
    # From table and text
    V_RF = 1.0  # V, peak amplitude of fundamental
    f_fundamental = 915e6  # Hz
    Cc = 1e-12  # F (1pF)
    RL = 8e3  # Ohm (8kΩ)
    CL = 5e-12  # F (5pF)
    
    # From problem description
    R0_parasitic = 50.0  # Ohm
    f0_parasitic = 915e6  # Hz
    C_parasitic = 2e-15  # F (2fF)
    
    harmonics = [1, 3, 5, 7]
    output_peak_voltages = []

    # Total load capacitance (CL and C_parasitic are in parallel)
    C_load_total = CL + C_parasitic

    # --- 2. Frequency-by-Frequency Analysis ---
    for n in harmonics:
        # --- 3. Calculate Input Harmonic Amplitude ---
        # Voltage drops by 10% (multiplied by 0.9) for each higher harmonic order
        # from the previous one in the list (1->3, 3->5, 5->7)
        V_in_peak = V_RF * (0.9)**((n-1)/2)

        # --- 4.a Calculate frequency ---
        f_n = n * f_fundamental
        w_n = 2 * np.pi * f_n

        # --- 4.b Calculate frequency-dependent resistance ---
        # R_parasitic(f) = R0 * (f/f0)^2. Since f_fundamental = f0, this is R0*n^2
        R_p_n = R0_parasitic * (n**2)

        # --- 4.c Calculate complex impedances ---
        # Impedance of coupling capacitor Cc
        Z_Cc_n = 1 / (1j * w_n * Cc)

        # --- 4.d Calculate combined load impedance ---
        # Admittance of the parallel load (1/RL + jw(CL+Cp))
        Y_L_n = 1/RL + 1j * w_n * C_load_total
        Z_L_n = 1 / Y_L_n
        
        # --- 4.e Use voltage divider to find output voltage ---
        # V_out = V_in * Z_L / (Z_series + Z_L)
        Z_total_series = R_p_n + Z_Cc_n + Z_L_n
        V_out_complex = V_in_peak * Z_L_n / Z_total_series

        # --- 4.f Get the magnitude (peak amplitude) ---
        V_out_peak = np.abs(V_out_complex)
        output_peak_voltages.append(V_out_peak)

    # --- 5. Calculate Total RMS Voltage ---
    v_l1, v_l3, v_l5, v_l7 = output_peak_voltages
    
    # Sum of the squares of the peak voltages
    sum_of_squares = v_l1**2 + v_l3**2 + v_l5**2 + v_l7**2
    
    # Total RMS voltage is sqrt( sum(V_peak_n^2 / 2) )
    V_total_rms = np.sqrt(sum_of_squares / 2)

    # --- 6. Final Output ---
    print(f"Calculated peak voltages at the load for each harmonic:")
    print(f"|V_L_1| = {v_l1:.4f} V")
    print(f"|V_L_3| = {v_l3:.4f} V")
    print(f"|V_L_5| = {v_l5:.4f} V")
    print(f"|V_L_7| = {v_l7:.4f} V\n")

    print(f"The final equation for the total RMS voltage is:")
    print(f"V_total_rms = sqrt( (|V_L_1|^2 + |V_L_3|^2 + |V_L_5|^2 + |V_L_7|^2) / 2 )")
    print(f"V_total_rms = sqrt( ( {v_l1:.4f}^2 + {v_l3:.4f}^2 + {v_l5:.4f}^2 + {v_l7:.4f}^2 ) / 2 )")
    print(f"V_total_rms = sqrt( ( {v_l1**2:.4f} + {v_l3**2:.4f} + {v_l5**2:.4f} + {v_l7**2:.4f} ) / 2 )")
    print(f"V_total_rms = sqrt( {sum_of_squares:.4f} / 2 ) = {V_total_rms:.4f} V")
    
    # Final answer for the platform
    # To return only the number as requested
    # print(f"<<<{V_total_rms:.4f}>>>")

if __name__ == '__main__':
    calculate_total_voltage()
<<<0.4027>>>