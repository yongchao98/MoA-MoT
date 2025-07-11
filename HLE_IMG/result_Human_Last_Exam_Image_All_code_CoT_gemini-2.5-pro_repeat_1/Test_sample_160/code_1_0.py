import numpy as np

def solve_circuit_voltage():
    """
    Calculates the total voltage across the load capacitor C_L,
    considering harmonic distortion and parasitic effects.
    """
    # --- Step 1: Define Constants and Initial Harmonic Amplitudes ---
    f_fundamental = 915e6  # 915 MHz
    V_rf = 1.0             # 1V peak for the fundamental
    R_L = 8000.0           # 8 kOhm
    C_parasitic = 2e-15    # 2 fF
    R0_parasitic = 50.0    # 50 Ohm for parasitic R calculation
    f0_parasitic = 915e6   # 915 MHz for parasitic R calculation

    harmonics = [1, 3, 5, 7]
    V_s_peaks = {
        1: V_rf,
        3: V_rf * 0.9,
        5: V_rf * 0.9**2,
        7: V_rf * 0.9**3
    }

    # --- Step 2: Model the Rectifier's Input Impedance ---
    # A common approximation for R_in of a cross-coupled rectifier is R_L / 2
    Z_in_rectifier = R_L / 2.0

    print("### Calculation Plan ###")
    print("1. Calculate attenuation for each harmonic due to parasitics.")
    print("2. Determine the peak voltage of each harmonic at the rectifier's input.")
    print(f"3. Combine harmonic peaks using RMS sum: V_total = sqrt(V_1^2 + V_3^2 + V_5^2 + V_7^2).")
    print(f"4. Assumed rectifier input resistance: R_in = R_L/2 = {Z_in_rectifier:.0f} Ohms\n")

    V_rect_in_peak_sq = []
    
    # --- Step 3 & 4: Calculate for each harmonic and prepare for summation ---
    print("### Detailed Harmonic Calculations ###")
    for n in harmonics:
        f_n = n * f_fundamental
        w_n = 2 * np.pi * f_n
        V_s_n = V_s_peaks[n]

        # Calculate frequency-dependent parasitics
        R_p_n = R0_parasitic * (f_n / f0_parasitic)**2
        Z_Cp_n = 1 / (1j * w_n * C_parasitic)

        # Calculate the effective load impedance seen by the source
        # (rectifier input impedance in parallel with parasitic capacitance)
        Z_shunt_n = (Z_in_rectifier * Z_Cp_n) / (Z_in_rectifier + Z_Cp_n)

        # Calculate the voltage divider attenuation factor
        # This is | V_in_rectifier / V_source |
        attenuation_factor = np.abs(Z_shunt_n / (R_p_n + Z_shunt_n))

        # Calculate the peak voltage of the harmonic at the rectifier input
        V_rect_in_n = V_s_n * attenuation_factor
        V_rect_in_peak_sq.append(V_rect_in_n**2)

        print(f"--- Harmonic n={n} ---")
        print(f"Source Voltage (V_S_{n}): {V_s_n:.3f} V")
        print(f"Parasitic Resistance (R_p_{n}): {R_p_n:.1f} \u03A9")
        print(f"Voltage at Rectifier Input (V_rect_{n}):")
        print(f"  V_rect_{n} = {V_s_n:.3f} * | ({Z_shunt_n.real:.1f} + {Z_shunt_n.imag:.1f}j) / ({R_p_n:.1f} + ({Z_shunt_n.real:.1f} + {Z_shunt_n.imag:.1f}j)) |")
        print(f"  V_rect_{n} = {V_s_n:.3f} * {attenuation_factor:.4f} = {V_rect_in_n:.4f} V\n")

    # --- Step 5: Calculate Final Total Voltage ---
    total_voltage = np.sqrt(np.sum(V_rect_in_peak_sq))
    
    print("### Final Voltage Calculation ###")
    equation_str = "V_total = sqrt("
    for i, val in enumerate(V_rect_in_peak_sq):
        v_peak = np.sqrt(val)
        equation_str += f"{v_peak:.4f}\u00b2"
        if i < len(V_rect_in_peak_sq) - 1:
            equation_str += " + "
    equation_str += ")"
    print(f"Equation: {equation_str}")
    
    equation_str_vals = "V_total = sqrt("
    for i, val in enumerate(V_rect_in_peak_sq):
        equation_str_vals += f"{val:.4f}"
        if i < len(V_rect_in_peak_sq) - 1:
            equation_str_vals += " + "
    equation_str_vals += ")"
    print(f"Values:   {equation_str_vals}")
    print(f"Result:   V_total = sqrt({np.sum(V_rect_in_peak_sq):.4f}) = {total_voltage:.4f} V")
    
    return total_voltage

if __name__ == '__main__':
    final_voltage = solve_circuit_voltage()
    print(f"\nThe calculated total voltage across C_L is: {final_voltage:.4f} V")
    # The final answer is wrapped in <<<>>> as requested.
    print(f"\n<<<{final_voltage:.4f}>>>")
