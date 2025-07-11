import math

def calculate_total_voltage():
    """
    Calculates the total voltage across the load capacitor CL, considering
    harmonic distortions and parasitic effects.
    """
    # --- Constants from the problem description and circuit diagram ---
    V_RF = 1.0  # V, Fundamental peak voltage
    f_fundamental = 915e6  # Hz, Fundamental frequency
    R_L = 8e3  # Ohm
    C_L = 5e-12  # F (5 pF)
    
    # Parasitic and frequency-dependent parameters
    R0 = 50.0  # Ohm
    f0 = 915e6  # Hz, Reference frequency for parasitic resistance
    C_parasitic = 2e-15  # F (2 fF)
    
    # Harmonics to consider
    harmonics = [1, 3, 5, 7]
    
    # Total effective capacitance at the load
    C_total = C_L + C_parasitic
    
    total_voltage_dc = 0.0
    V_in_amp = V_RF
    harmonic_voltages = []

    print("--- Calculating DC Voltage Contribution from Each Harmonic ---")
    
    # Calculate for each harmonic
    for n in harmonics:
        if n > 1:
            # Voltage drops by 10% for each higher harmonic
            V_in_amp *= 0.9
            
        # 1. Calculate frequency for the current harmonic
        f_n = n * f_fundamental
        
        # 2. Calculate frequency-dependent parasitic resistance
        R_p = R0 * (f_n / f0)**2
        
        # 3. Calculate capacitive reactance of the total load capacitance
        # Xc = 1 / (2 * pi * f * C)
        Xc_total = 1 / (2 * math.pi * f_n * C_total)
        
        # 4. Calculate the attenuated voltage amplitude using the R-C divider model
        # V_attenuated = V_in * |Zc| / sqrt(R^2 + Zc^2)
        V_attenuated = V_in_amp * Xc_total / math.sqrt(R_p**2 + Xc_total**2)
        
        # 5. Calculate DC voltage contribution assuming an ideal voltage doubler (gain=2)
        V_dc_n = 2 * V_attenuated
        harmonic_voltages.append(V_dc_n)
        
        # Print the detailed calculation for the current harmonic
        print(f"\nFor Harmonic n={n}:")
        print(f"  - Input Voltage (V_in_{n}): {V_in_amp:.3f} V")
        print(f"  - Frequency (f_{n}): {f_n/1e6:.0f} MHz")
        print(f"  - Parasitic Resistance (R_p_{n}): {R_p:.2f} \u03A9")
        print(f"  - Total Capacitive Reactance (X_c_{n}): {Xc_total:.2f} \u03A9")
        print(f"  - Attenuated Voltage: {V_attenuated:.4f} V")
        print(f"  - DC Contribution (V_dc_{n}): 2 * {V_attenuated:.4f} V = {V_dc_n:.4f} V")

    # Sum the DC contributions
    total_voltage_dc = sum(harmonic_voltages)
    
    # Print the final equation
    print("\n--- Final Calculation ---")
    equation_parts = [f"{v:.4f} V" for v in harmonic_voltages]
    equation = " + ".join(equation_parts)
    print(f"Total DC Voltage = {equation} = {total_voltage_dc:.4f} V")

# Execute the calculation
calculate_total_voltage()

<<<1.2003>>>