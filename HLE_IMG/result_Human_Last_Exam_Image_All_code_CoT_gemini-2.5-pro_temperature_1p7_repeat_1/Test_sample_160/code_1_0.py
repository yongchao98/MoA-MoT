import numpy as np

def solve_circuit_voltage():
    """
    Calculates the total voltage across the load capacitor CL, considering
    harmonic distortions and parasitic effects.
    """

    # --- Step 1: Define Constants from the problem statement and image ---
    V_RF = 1.0  # Peak voltage of the fundamental frequency in Volts
    f0 = 915e6  # Fundamental frequency in Hz (from table omega = 2*pi*915MHz)
    R_L = 8e3     # Load resistance in Ohms
    C_L = 5e-12     # Load capacitance in Farads
    R0 = 50.0    # Base parasitic resistance in Ohms for the formula
    f_ref = 915e6 # Reference frequency for parasitic resistance in Hz
    C_parasitic = 2e-15 # Parasitic capacitance in Farads

    # --- Step 2: Characterize the input signal with harmonics ---
    harmonics = [1, 3, 5, 7]
    V_in_peak = {}
    V_in_peak[1] = V_RF
    # The voltage drops by 10% for each higher harmonic relative to the previous one
    for i in range(1, len(harmonics)):
        prev_harmonic = harmonics[i-1]
        current_harmonic = harmonics[i]
        V_in_peak[current_harmonic] = V_in_peak[prev_harmonic] * 0.9

    # --- Step 3: Model parasitic and load effects ---
    # The total capacitance at the load is CL in parallel with C_parasitic
    C_total_load = C_L + C_parasitic

    V_out_phasors = {}
    V_out_amplitudes = {}
    
    # --- Step 4: Calculate the transfer function for each harmonic ---
    for k in harmonics:
        f_k = k * f0
        w_k = 2 * np.pi * f_k

        # Frequency-dependent parasitic resistance
        R_parasitic_k = R0 * (f_k / f_ref)**2

        # Impedance of the total capacitance at the load
        Z_C_k = 1 / (1j * w_k * C_total_load)

        # Impedance of the parallel RL and C_total_load combination
        Z_load_k = (R_L * Z_C_k) / (R_L + Z_C_k)

        # Calculate the output voltage phasor using the voltage divider rule
        V_out_phasor = V_in_peak[k] * (Z_load_k / (R_parasitic_k + Z_load_k))
        V_out_phasors[k] = V_out_phasor
        V_out_amplitudes[k] = np.abs(V_out_phasor)
        
    # --- Step 5: Determine the peak voltage of the composite output signal ---
    # We create the time-domain signal by summing the contribution from each harmonic
    # and find its peak value numerically.
    T0 = 1/f0  # Period of the fundamental frequency
    # High-resolution time vector to capture peaks accurately
    t = np.linspace(0, T0, 5000) 
    v_out_t = np.zeros_like(t)

    for k in harmonics:
        w_k = 2 * np.pi * k * f0
        phasor = V_out_phasors[k]
        amplitude = np.abs(phasor)
        phase = np.angle(phasor)
        # We use cosine since phasors are defined relative to cos(wt)
        v_out_t += amplitude * np.cos(w_k * t + phase)

    # The total voltage is the peak of this composite waveform
    V_total_peak = np.max(v_out_t)

    # --- Step 6: Present the Final Answer ---
    print("The total voltage across the load capacitor is the peak of the composite output waveform.")
    print("The output voltage waveform is a sum of harmonics:")
    print("v_out(t) = V_1*cos(w1*t + p1) + V_3*cos(w3*t + p3) + V_5*cos(w5*t + p5) + V_7*cos(w7*t + p7)")
    print("\nThe peak amplitudes (V_k) of the output voltage harmonics are:")
    for k in sorted(V_out_amplitudes.keys()):
        print(f"V_{k} = {V_out_amplitudes[k]:.4f} V")

    print("\nAfter combining these harmonic components, the final equation is:")
    equation_parts = []
    for k in sorted(V_out_amplitudes.keys()):
        equation_parts.append(f"{V_out_amplitudes[k]:.4f} * cos({k}*w*t + phi_{k})")
    final_equation = "v_out(t) = " + " + ".join(equation_parts)
    print(final_equation)
    
    print(f"\nThe total voltage across C_L (peak of v_out(t)) is: {V_total_peak:.4f} V")
    
    return V_total_peak

# Execute the function to get the result
solve_circuit_voltage()
<<<0.5750>>>