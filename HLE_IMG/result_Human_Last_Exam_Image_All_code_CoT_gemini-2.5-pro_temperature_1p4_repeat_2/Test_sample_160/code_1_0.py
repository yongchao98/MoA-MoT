import math

def calculate_output_voltage():
    """
    Calculates the total DC voltage across the load capacitor CL, considering
    harmonic distortion and parasitic effects.
    """
    # --- Parameters from the problem description and image ---
    # Source
    V_RF_peak = 1.0  # V, peak voltage of the fundamental
    f1 = 915e6       # Hz, fundamental frequency

    # Harmonics
    harmonic_orders = [1, 3, 5, 7]
    harmonic_attenuation = 0.90  # 10% drop relative to the previous harmonic

    # Circuit components
    Cc = 1e-12       # F
    RL = 8e3         # Ohms

    # Parasitics
    R0 = 50.0        # Ohms, base parasitic resistance
    f0 = 915e6       # Hz, reference frequency for parasitic resistance

    # --- Model Assumptions ---
    # Effective input capacitance of the rectifier (two Cc in series)
    C_in_eff = Cc / 2.0
    # Ideal voltage gain of the voltage doubler rectifier
    G_v = 2.0

    # --- Calculations ---
    V_source_harmonics = {}

    # Calculate source voltage for each harmonic
    V_current = V_RF_peak
    for i, n in enumerate(harmonic_orders):
        if i == 0:
            V_source_harmonics[n] = V_RF_peak
        else:
            # Voltage is based on the previous harmonic in the sequence
            prev_n = harmonic_orders[i-1]
            V_source_harmonics[n] = V_source_harmonics[prev_n] * harmonic_attenuation

    print("Calculating DC voltage contribution from each harmonic...\n")

    total_V_DC = 0.0
    equation_parts = []

    for n in harmonic_orders:
        f_n = n * f1
        V_n = V_source_harmonics[n]

        # a. Calculate frequency-dependent parasitic resistance
        R_p_n = R0 * (f_n / f0)**2

        # b. Calculate the magnitude of the rectifier's input impedance (capacitive)
        Xc_n = 1 / (2 * math.pi * f_n * C_in_eff)

        # c. Calculate the peak voltage at the rectifier input (voltage divider)
        V_rect_in_n = V_n * Xc_n / math.sqrt(R_p_n**2 + Xc_n**2)

        # d. Calculate the unloaded DC output voltage from this harmonic
        V_DC_unloaded_n = G_v * V_rect_in_n

        # e. Calculate the rectifier's approximate output resistance
        R_out_rect_n = 1 / (f_n * Cc)

        # f. Calculate the loaded DC output voltage (voltage divider with RL)
        V_DC_loaded_n = V_DC_unloaded_n * RL / (R_out_rect_n + RL)
        
        total_V_DC += V_DC_loaded_n
        equation_parts.append(f"{V_DC_loaded_n:.4f}")
        
        print(f"--- Harmonic n={n} ---")
        print(f"  Loaded DC Voltage Contribution (V_DC_{n}): {V_DC_loaded_n:.4f} V")
        print("-" * 22)

    # --- Final Result ---
    final_equation = " + ".join(equation_parts)
    print("\nThe total voltage across C_L is the sum of the DC contributions from each harmonic.")
    print(f"Total V_DC = V_DC_1 + V_DC_3 + V_DC_5 + V_DC_7")
    print(f"Total V_DC = {final_equation} = {total_V_DC:.4f} V")

# Execute the calculation
calculate_output_voltage()
<<<2.2875>>>