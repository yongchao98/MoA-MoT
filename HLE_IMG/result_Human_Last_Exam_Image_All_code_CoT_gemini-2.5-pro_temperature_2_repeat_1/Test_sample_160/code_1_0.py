import numpy as np

def calculate_total_voltage():
    """
    Calculates the total voltage across the load capacitor C_L, considering
    harmonic distortions, parasitic effects, and frequency-dependent losses.
    """
    # 1. Define Constants and Inputs from the problem description
    V_RF = 1.0  # V, peak amplitude of the fundamental frequency
    f_fund = 915e6  # Hz, fundamental frequency from w = 2*pi*915MHz
    C_c = 1e-12  # F, coupling capacitance
    R_L = 8e3  # Ohms, load resistance
    C_L = 5e-12  # F, load capacitance
    R0 = 50.0  # Ohms, base parasitic resistance
    f0_ref = 915e6  # Hz, reference frequency for parasitic resistance
    C_parasitic = 2e-15  # F, parasitic capacitance

    # Harmonic orders to consider
    harmonic_orders = [1, 3, 5, 7]

    # 2. Model Harmonic Amplitudes
    # Assume "significant" 3rd harmonic distortion is 10% of the fundamental.
    # Subsequent harmonics drop by 10% (i.e., are 90% of the previous).
    V_in_amplitudes = {}
    V_in_amplitudes[1] = V_RF
    # HD3 is 10% of fundamental
    V_in_amplitudes[3] = 0.1 * V_in_amplitudes[1]
    # HD5 is 90% of HD3
    V_in_amplitudes[5] = 0.9 * V_in_amplitudes[3]
    # HD7 is 90% of HD5
    V_in_amplitudes[7] = 0.9 * V_in_amplitudes[5]

    V_out_amplitudes = []
    
    # 3. AC Circuit Analysis for Each Harmonic
    for n in harmonic_orders:
        f = n * f_fund
        w = 2 * np.pi * f
        V_in = V_in_amplitudes[n]

        # 4. Calculate Component Impedances
        # Parasitic resistance (frequency-dependent)
        R_p = R0 * (f / f0_ref)**2
        
        # Series impedance of coupling capacitor
        Z_Cc = 1 / (1j * w * C_c)
        
        # Total load capacitance
        C_L_total = C_L + C_parasitic
        
        # Impedance of the total load capacitance
        Z_CL_total = 1 / (1j * w * C_L_total)
        
        # Equivalent impedance of the parallel R_L and C_L_total
        Z_load = (R_L * Z_CL_total) / (R_L + Z_CL_total)
        
        # Total series impedance
        Z_series = R_p + Z_Cc

        # 5. Calculate Output Voltage for the harmonic
        V_out_complex = V_in * Z_load / (Z_series + Z_load)
        V_out_peak = np.abs(V_out_complex)
        V_out_amplitudes.append(V_out_peak)

    # 6. Calculate Total Voltage by summing the peak amplitudes
    total_voltage = sum(V_out_amplitudes)

    # 7. Final Output
    # The problem asks to output the final equation with each number.
    print("The total voltage across C_L is the sum of the peak voltage contributions from each harmonic.")
    equation_str = "V_total = V_out(h1) + V_out(h3) + V_out(h5) + V_out(h7)"
    print(equation_str)
    
    value_str = "V_total = {:.4f} V + {:.4f} V + {:.4f} V + {:.4f} V".format(
        V_out_amplitudes[0], V_out_amplitudes[1], V_out_amplitudes[2], V_out_amplitudes[3]
    )
    print(value_str)
    
    result_str = "V_total = {:.4f} V".format(total_voltage)
    print(result_str)
    
    # Returning the final numerical value for the grading system
    return total_voltage

# Execute the calculation and print the result
final_voltage = calculate_total_voltage()
# The final answer format for the platform
print(f"<<<{final_voltage:.4f}>>>")