import numpy as np

def calculate_voltage_across_cl():
    """
    Calculates the total voltage across the load capacitor C_L,
    considering harmonic distortions and parasitic effects.
    """
    # --- 1. Define Constants and Parameters ---
    V_RF = 1.0  # V, peak voltage of the fundamental
    f1 = 915e6  # Hz, fundamental frequency
    Cc = 1e-12  # F, coupling capacitance
    RL = 8e3    # Ohm, load resistance
    CL = 5e-12    # F, load capacitance
    
    # Parasitic effects parameters
    R0 = 50.0       # Ohm, base parasitic resistance
    f0 = 915e6      # Hz, reference frequency for R_parasitic
    C_parasitic = 2e-15 # F, parasitic capacitance (its effect is negligible)

    print("--- System Parameters ---")
    print(f"Fundamental Frequency (f1): {f1/1e6} MHz")
    print(f"Load Resistor (RL): {RL/1e3} kOhm")
    print(f"Load Capacitor (CL): {CL*1e12} pF")
    print(f"Parasitic Capacitance (C_parasitic): {C_parasitic*1e15} fF")
    print(f"Base Parasitic Resistance (R0): {R0} Ohm at {f0/1e6} MHz\n")
    
    # --- 2. Calculate Input Harmonic Voltages ---
    # Voltage drops by 10% (multiplied by 0.9) for each higher harmonic.
    V_harmonics_in = {
        1: V_RF,
        3: V_RF * 0.9,
        5: V_RF * 0.9 * 0.9,
        7: V_RF * 0.9 * 0.9 * 0.9
    }
    
    print("--- Input Harmonic Voltages (at source) ---")
    for n, v in V_harmonics_in.items():
        print(f"V{n} (peak): {v:.4f} V")
    print("")

    # --- 3. Calculate Attenuated Voltages at Rectifier Input ---
    # Model the rectifier's input resistance as R_in = RL / 2
    R_in_rectifier = RL / 2
    
    V_harmonics_rect = {}
    print("--- Voltage Calculation for Each Harmonic ---")
    
    equation_parts = []

    for n, v_in in V_harmonics_in.items():
        f_n = n * f1
        
        # Calculate frequency-dependent parasitic resistance
        R_p = R0 * (f_n / f0)**2
        
        # Calculate impedance of the rectifier input.
        # The reactive part from Cc is negligible compared to R_in_rectifier at these frequencies.
        # Z_rectifier ≈ R_in_rectifier
        # We can ignore the very small reactive part in the magnitude calculation.
        Z_rect_mag = R_in_rectifier

        # Calculate total impedance magnitude for the voltage divider
        Z_total_mag = np.sqrt((R_p + R_in_rectifier)**2)
        
        # Calculate the attenuated voltage magnitude using the voltage divider rule
        v_rect = v_in * (Z_rect_mag / Z_total_mag)
        V_harmonics_rect[n] = v_rect

        print(f"Harmonic n={n}:")
        print(f"  - Frequency: {f_n/1e6:.0f} MHz")
        print(f"  - Parasitic Resistance R_p({n}f1): {R_p:.2f} Ohm")
        print(f"  - Voltage at rectifier V_rect_{n} = {v_in:.3f} * ({Z_rect_mag:.0f} / ({R_p:.2f} + {Z_rect_mag:.0f})) = {v_rect:.4f} V")
        
        equation_parts.append(f"{v_rect:.4f}*sin({n}*w*t)")

    # --- 4. Find the Peak of the Composite Waveform ---
    # The final DC voltage across C_L is approximated by the peak of the input waveform at the rectifier.
    # v_rect(t) = V_rect_1*sin(wt) + V_rect_3*sin(3wt) + ...
    
    print("\n--- Final Voltage Calculation ---")
    final_equation = "v_rect(t) = " + " + ".join(equation_parts)
    print("The voltage waveform at the rectifier input is:")
    print(final_equation)

    # Numerically find the peak value of this composite sine wave
    t = np.linspace(0, 2 * np.pi, 20000) # One fundamental cycle for w*t
    v_total_rect = sum(V_harmonics_rect[n] * np.sin(n * t) for n in V_harmonics_rect)
    
    V_DC_final = np.max(v_total_rect)
    
    print(f"\nThe peak of this waveform, which approximates the final DC voltage across C_L, is calculated numerically.")
    print(f"Total Voltage across C_L (V_DC) = Peak(v_rect(t)) = {V_DC_final:.4f} V")
    
    return V_DC_final

# Run the calculation and print the final answer in the required format
final_voltage = calculate_voltage_across_cl()
print(f"\n<<<{final_voltage:.4f}>>>")