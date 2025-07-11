import numpy as np

def calculate_total_voltage():
    """
    Calculates the total RMS voltage across the load capacitor C_L,
    considering harmonic distortions and parasitic effects.
    """
    # --- 1. Define Constants from the problem description and image ---
    V_RF_peak = 1.0         # V (Amplitude of the fundamental input signal)
    f_fundamental = 915e6   # Hz (Fundamental frequency)
    R_L = 8e3               # Ohms (Load resistance)
    C_L = 5e-12             # Farads (Load capacitance)

    # Parasitic component parameters
    R0_parasitic = 50.0     # Ohms (Base parasitic resistance)
    f0_parasitic = 915e6    # Hz (Reference frequency for parasitic resistance)
    C_parasitic = 2e-15     # Farads (Parasitic capacitance)

    # Harmonic information
    harmonics = [1, 3, 5, 7]
    distortion_factor = 0.9 # Voltage drops by 10% for each higher harmonic

    # --- 2. Perform Calculations ---
    output_peak_voltages_sq = []
    output_peak_voltages = []
    V_in_peak = V_RF_peak

    print("Calculating total voltage across C_L...")
    print("-" * 70)
    print(f"Fundamental Frequency (f1): {f_fundamental/1e6:.0f} MHz")
    print(f"Input Voltage Amplitude (V_in1): {V_RF_peak:.2f} V")
    print(f"Load Resistance (R_L): {R_L/1e3:.1f} kOhm")
    print(f"Load Capacitance (C_L): {C_L*1e12:.1f} pF")
    print(f"Parasitic Capacitance (C_parasitic): {C_parasitic*1e15:.1f} fF")
    print("-" * 70)

    # The total capacitance at the load is C_L in parallel with C_parasitic
    C_total = C_L + C_parasitic

    for n in harmonics:
        # For harmonics n > 1, the input voltage is reduced
        if n > 1:
            V_in_peak *= distortion_factor

        # Calculate frequency for the current harmonic
        f_n = n * f_fundamental
        omega_n = 2 * np.pi * f_n

        # Calculate frequency-dependent parasitic resistance
        R_p_n = R0_parasitic * (f_n / f0_parasitic)**2

        # Calculate impedance of the total load capacitance
        # Z = 1 / (j*omega*C)
        Z_C_total_n = 1 / (1j * omega_n * C_total)

        # Calculate the parallel impedance of R_L and the total capacitance
        # Z_parallel = (Z1 * Z2) / (Z1 + Z2)
        Z_L_parallel_n = (R_L * Z_C_total_n) / (R_L + Z_C_total_n)

        # Calculate the output voltage amplitude for this harmonic using the voltage divider rule
        # V_out = V_in * Z_load / (Z_series + Z_load)
        V_out_peak_n = V_in_peak * abs(Z_L_parallel_n / (R_p_n + Z_L_parallel_n))

        output_peak_voltages.append(V_out_peak_n)
        output_peak_voltages_sq.append(V_out_peak_n**2)

        print(f"For Harmonic n = {n}:")
        print(f"  Input Voltage Peak: {V_in_peak:.4f} V")
        print(f"  Frequency: {f_n/1e6:.2f} MHz")
        print(f"  Parasitic Resistance: {R_p_n:.4f} Ohms")
        print(f"  Output Voltage Peak: {V_out_peak_n:.4f} V\n")

    # --- 3. Calculate Total RMS Voltage ---
    # V_total_rms = sqrt( V_rms1^2 + V_rms3^2 + ... )
    # Since V_rms_n = V_peak_n / sqrt(2), then V_total_rms^2 = sum(V_peak_n^2 / 2)
    # So, V_total_rms = sqrt( sum(V_peak_n^2) ) / sqrt(2)
    total_rms_voltage = np.sqrt(sum(output_peak_voltages_sq)) / np.sqrt(2)

    # --- 4. Final Output ---
    print("-" * 70)
    print("The total RMS voltage is calculated from the peak voltages of each harmonic:")
    print("V_total_rms = (1/sqrt(2)) * sqrt(V_peak1^2 + V_peak3^2 + V_peak5^2 + V_peak7^2)")
    print("\nFinal Equation:")
    
    v_p1, v_p3, v_p5, v_p7 = output_peak_voltages
    equation_str = (
        f"V_total_rms = (1/sqrt(2)) * sqrt({v_p1:.4f}^2 + {v_p3:.4f}^2 + "
        f"{v_p5:.4f}^2 + {v_p7:.4f}^2)"
    )
    print(equation_str)
    
    print(f"\nResult:")
    print(f"V_total_rms = {total_rms_voltage:.4f} V")

    return total_rms_voltage

if __name__ == '__main__':
    result = calculate_total_voltage()
    print(f"<<<{result:.4f}>>>")