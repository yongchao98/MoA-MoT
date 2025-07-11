import math
import cmath

def solve_circuit_voltage():
    """
    Calculates the total voltage across the load in the given circuit,
    considering harmonic distortions and parasitic effects.
    """
    # 1. Define constants from the problem description and image
    V_RF = 1.0         # V, peak amplitude of the fundamental
    f1 = 915e6         # Hz, fundamental frequency
    R0 = 50.0          # Ohms, base parasitic resistance
    f0 = 915e6         # Hz, reference frequency for parasitic resistance
    C_parasitic = 2e-15 # F, (2 fF)
    RL = 8e3           # Ohms, (8 kOhm)
    CL = 5e-12         # F, (5 pF)

    # Combine load and parasitic capacitance
    C_total = CL + C_parasitic

    # List of harmonics to consider
    harmonics = [1, 3, 5, 7]
    
    # Dictionary to store results for each harmonic
    results = {}

    print("Calculating total voltage across the load...")
    print("-" * 40)

    # 2. Loop through each harmonic to calculate its contribution
    V_in_peak = V_RF
    for n in harmonics:
        # Calculate frequency-dependent values
        fn = n * f1
        omega_n = 2 * math.pi * fn
        
        # Update input voltage for harmonics > 1
        if n > 1:
            # Voltage drops by 10% for each higher harmonic
            V_in_peak *= 0.9

        # Calculate series parasitic resistance for the current harmonic
        R_s_n = R0 * (fn / f0)**2

        # Calculate complex load impedance Z_L = R_L || (C_L + C_parasitic)
        # Z_L = 1 / (1/R_L + j*omega*C_total)
        Z_L_n = 1 / (1/RL + 1j * omega_n * C_total)

        # Use the voltage divider formula to find the output voltage across the load
        # V_out = V_in * Z_L / (R_s + Z_L)
        V_out_complex = V_in_peak * (Z_L_n / (R_s_n + Z_L_n))
        V_out_peak = abs(V_out_complex)

        # Store results
        results[n] = {
            'V_in_peak': V_in_peak,
            'R_s': R_s_n,
            'V_out_peak': V_out_peak
        }

        # Print intermediate results for clarity
        print(f"Harmonic n={n}:")
        print(f"  - Input Voltage (peak): {V_in_peak:.3f} V")
        print(f"  - Parasitic Resistance: {R_s_n:.1f} Ohm")
        print(f"  - Output Voltage (peak): {V_out_peak:.4f} V\n")

    # 3. Calculate the total voltage
    # The final DC voltage is assumed to be the total RMS voltage of all harmonic
    # components across the load, which is sqrt( sum(V_peak^2) / 2 ).
    sum_V_out_sq = sum(res['V_out_peak']**2 for res in results.values())
    V_total = math.sqrt(sum_V_out_sq / 2)

    # 4. Print the final calculation and result
    print("-" * 40)
    print("The final DC voltage is the total RMS voltage of all harmonic components across the load.")
    
    # Build and print the equation string
    v_out_sq_terms = " + ".join(f"{results[n]['V_out_peak']:.4f}^2" for n in harmonics)
    print(f"V_total = sqrt( ({v_out_sq_terms}) / 2 )")
    print(f"V_total = {V_total:.4f} V")

# Execute the function
solve_circuit_voltage()
<<<0.4028>>>