import numpy as np
import cmath

def solve_circuit_voltage():
    """
    Calculates the total voltage across the load capacitor CL in a rectifier circuit,
    considering harmonics and parasitic effects.
    """
    # --- Parameters from the image and problem description ---
    V_RF = 1.0  # V, amplitude of the fundamental
    f1 = 915e6  # Hz, fundamental frequency
    RL = 8e3    # Ohms
    CL = 5e-12  # F
    R0 = 50.0   # Ohms, for parasitic resistance calculation
    f0 = 915e6  # Hz, reference frequency for parasitic resistance
    C_parasitic = 2e-15  # F (2 fF)
    harmonic_drop_factor = 0.9

    # The total capacitance at the load is the sum of the load and parasitic capacitances.
    C_total = CL + C_parasitic

    # --- Step 1: Define input harmonic amplitudes ---
    harmonics_in = {
        1: V_RF,
        3: V_RF * harmonic_drop_factor,
        5: V_RF * harmonic_drop_factor**2,
        7: V_RF * harmonic_drop_factor**3
    }

    print("--- Input Signal and System Parameters ---")
    print(f"Fundamental Frequency f1 = {f1/1e6:.0f} MHz")
    print(f"Load Resistor RL = {RL/1e3:.1f} kOhm")
    print(f"Total Load Capacitance C_total = {C_total*1e12:.4f} pF")
    print("Initial Harmonic Amplitudes (V_in):")
    for n, V_in in harmonics_in.items():
        print(f"  V_{n} = {V_in:.4f} V")
    print("-" * 40)

    # --- Step 2: Analyze the circuit for each harmonic ---
    output_phasors = {}
    print("\n--- Frequency-Dependent Filtering Analysis ---")

    for n, V_in in harmonics_in.items():
        f_n = n * f1
        w_n = 2 * np.pi * f_n
        
        # Calculate frequency-dependent parasitic resistance
        Rs = R0 * (f_n / f0)**2
        
        # Calculate load impedance Zl = RL || (1/jwC_total)
        Zl = RL / (1 + 1j * w_n * RL * C_total)
        
        # Calculate complex transfer function H(f) = Zl / (Rs + Zl)
        H = Zl / (Rs + Zl)
        
        # Calculate the output phasor for the harmonic
        V_out_phasor = H * V_in
        output_phasors[n] = V_out_phasor

        print(f"\nProperties for Harmonic n={n} (f={f_n/1e6:.0f} MHz):")
        print(f"  Parasitic Resistance Rs = {Rs:.2f} Ohms")
        print(f"  Filtered Amplitude |V_out_{n}| = {abs(V_out_phasor):.4f} V")
        print(f"  Filtered Phase = {np.angle(V_out_phasor, deg=True):.2f} degrees")
    print("-" * 40)

    # --- Step 3: Find the peak of the output waveform ---
    T1 = 1 / f1
    num_points = 10000
    t = np.linspace(0, T1, num_points)
    w1 = 2 * np.pi * f1
    v_out_t = np.zeros_like(t)

    for n, V_phasor in output_phasors.items():
        v_out_t += abs(V_phasor) * np.cos(n * w1 * t + np.angle(V_phasor))
    
    V_peak_filtered = np.max(v_out_t)

    # --- Step 4: Calculate final voltage across CL ---
    # For a voltage doubler, V_DC ≈ 2 * V_peak
    V_CL = 2 * V_peak_filtered

    print("\n--- Final Voltage Calculation ---")
    print("The peak voltage of the filtered AC waveform (V_peak_filtered) is found by summing the harmonic components and finding the maximum.")
    
    # Print the equation for the final calculation
    print(f"\n1. Peak voltage of filtered signal: V_peak_filtered = {V_peak_filtered:.4f} V")
    print(f"2. The final voltage across C_L is calculated using the voltage doubler formula: V_CL = 2 * V_peak_filtered")
    print(f"   V_CL = 2 * {V_peak_filtered:.4f} V")
    print(f"   V_CL = {V_CL:.4f} V")
    
    # The final answer to be extracted
    # The format requires a clear final answer string.
    # print(f"\n<<<{V_CL:.4f}>>>")

solve_circuit_voltage()