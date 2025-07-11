import math

def calculate_total_voltage():
    """
    Calculates the total voltage across the load capacitor C_L, considering
    harmonic distortion and parasitic effects.
    """
    # 1. Define constants from the problem description and circuit diagram.
    V_RF = 1.0  # V, peak voltage of the fundamental frequency
    f1 = 915e6  # Hz, fundamental frequency
    R0 = 50.0   # Ohm, base parasitic resistance at f0
    f0 = 915e6  # Hz, reference frequency for parasitic resistance
    RL = 8000.0 # Ohm, load resistance
    # The 2fF parasitic capacitance has a negligible effect on the input impedance
    # compared to R_in and is omitted from the main calculation for simplicity.
    
    # 2. Make a standard approximation for the rectifier's input resistance.
    # For a voltage doubler, R_in is approximately RL / 2.
    R_in = RL / 2.0

    # 3. Define the harmonics to consider.
    harmonics = [1, 3, 5, 7]
    
    # 4. Calculate the source voltage for each harmonic.
    # The voltage drops by 10% for each higher harmonic.
    V_source = {}
    V_source[1] = V_RF
    V_source[3] = V_source[1] * 0.90
    V_source[5] = V_source[3] * 0.90
    V_source[7] = V_source[5] * 0.90
    
    # 5. Calculate the effective peak voltage at the rectifier input for each harmonic.
    V_eff = {}
    sum_of_squares = 0.0

    print("Step-by-Step Calculation of Effective Input Voltages:\n")
    for n in harmonics:
        # Parasitic resistance is frequency-dependent: R_p(f) = R0 * (f/f0)^2
        # Since f_n = n * f1 and f1 = f0, R_p_n = R0 * n^2
        R_p_n = R0 * (n**2)
        
        # Calculate the effective voltage using a voltage divider model: V_eff = V_source * R_in / (R_p + R_in)
        v_eff_n = V_source[n] * (R_in / (R_p_n + R_in))
        V_eff[n] = v_eff_n
        
        sum_of_squares += v_eff_n**2
        
        print(f"For Harmonic n={n}:")
        print(f"  - Source Voltage (V_{n}): {V_source[n]:.3f} V")
        print(f"  - Parasitic Resistance (R_p_{n}): {R_p_n:.1f} Ω")
        print(f"  - Effective Input Voltage (V_eff_{n}) = V_{n} * R_in / (R_p_{n} + R_in)")
        print(f"    V_eff_{n} = {V_source[n]:.3f} V * {R_in:.1f} Ω / ({R_p_n:.1f} Ω + {R_in:.1f} Ω) = {v_eff_n:.4f} V\n")

    # 6. Combine the effective voltages using an RMS sum to find the total effective peak voltage.
    V_peak_total = math.sqrt(sum_of_squares)
    
    # 7. Calculate the final DC voltage using the ideal voltage doubler formula (V_DC = 2 * V_peak).
    V_DC = 2 * V_peak_total

    # 8. Print the final calculation and result.
    print("--------------------------------------------------")
    print("Final Calculation for Total Voltage Across C_L:")
    print("--------------------------------------------------")
    print("1. Combine effective peak voltages from all harmonics:")
    print(f"V_peak_total = sqrt( V_eff_1² + V_eff_3² + V_eff_5² + V_eff_7² )")
    print(f"V_peak_total = sqrt( {V_eff[1]:.4f}² + {V_eff[3]:.4f}² + {V_eff[5]:.4f}² + {V_eff[7]:.4f}² )")
    print(f"V_peak_total = sqrt( {V_eff[1]**2:.4f} + {V_eff[3]**2:.4f} + {V_eff[5]**2:.4f} + {V_eff[7]**2:.4f} )")
    print(f"V_peak_total = sqrt( {sum_of_squares:.4f} ) = {V_peak_total:.4f} V\n")
    
    print("2. Calculate final DC voltage using the voltage doubler formula (V_DC = 2 * V_peak_total):")
    print(f"V_DC = 2 * {V_peak_total:.4f} V")
    print(f"V_DC = {V_DC:.4f} V\n")
    
    # Final answer formatted as requested
    print(f"<<<{V_DC:.4f}>>>")

calculate_total_voltage()