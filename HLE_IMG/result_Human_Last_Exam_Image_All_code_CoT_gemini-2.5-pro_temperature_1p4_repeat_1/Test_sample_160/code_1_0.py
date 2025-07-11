import math

def solve_circuit_voltage():
    """
    Calculates the total voltage across the load capacitor CL in the given rectifier circuit,
    considering harmonics and parasitic effects.
    """

    # --- Step 1: Define constants and parameters from the problem ---
    V_RF = 1.0  # V, amplitude of the fundamental input signal
    f1 = 915e6  # Hz, fundamental frequency from w = 2*pi*915MHz
    Cc = 1e-12  # F, coupling capacitance
    C_parasitic = 2e-15  # F, parasitic capacitance
    R0 = 50.0   # Ohm, base parasitic resistance
    f0 = 915e6  # Hz, reference frequency for parasitic resistance

    # Harmonics to consider in the input signal
    harmonics = [1, 3, 5, 7]

    # --- Step 2: Calculate harmonic amplitudes based on the given rule ---
    # "The voltage drops by 10% for each higher harmonic relative to the previous harmonic."
    # This means V_n = V_{prev_n} * (1 - 0.10)
    V_amps = {}
    V_amps[1] = V_RF
    # Note: The problem describes harmonics relative to the *previous* harmonic.
    # V3 is relative to V1, V5 relative to V3, V7 relative to V5.
    V_current = V_RF
    for i in range(1, len(harmonics)):
        V_current *= 0.9
        V_amps[harmonics[i]] = V_current

    total_voltage_across_CL = 0.0
    contributions = []

    print("--- Calculation Breakdown ---")

    # --- Step 3: Loop through each harmonic, analyze, and find its contribution ---
    for n in harmonics:
        V_in_amp = V_amps[n]
        fn = n * f1
        omega_n = 2 * math.pi * fn

        # Calculate frequency-dependent parasitic resistance
        R_para = R0 * (fn / f0)**2

        # Calculate complex impedances
        # Impedance of parasitic capacitance: Z = 1 / (j*omega*C)
        Z_cpara = 1 / (1j * omega_n * C_parasitic)

        # Input impedance of the rectifier (modeled as two series Cc)
        # Z = 2 / (j*omega*C)
        Z_rect_in = 2 / (1j * omega_n * Cc)

        # Calculate total load impedance seen by the source and series parasitic R
        # This is the rectifier input impedance in parallel with the parasitic capacitance.
        Z_load = (Z_rect_in * Z_cpara) / (Z_rect_in + Z_cpara)

        # Calculate the AC voltage amplitude at the rectifier's input using the voltage divider formula
        V_rect_in_complex = V_in_amp * Z_load / (R_para + Z_load)
        V_rect_in_amp = abs(V_rect_in_complex)

        # --- Step 4: Assume ideal rectification and sum contributions ---
        # V_DC_n = V_rect_in_amp
        V_dc_n = V_rect_in_amp
        total_voltage_across_CL += V_dc_n
        contributions.append(V_dc_n)

        print(f"Harmonic n = {n}:")
        print(f"  Input Voltage Amplitude = {V_in_amp:.3f} V")
        print(f"  Frequency = {fn/1e9:.3f} GHz")
        print(f"  Parasitic Resistance (R_para) = {R_para:.2f} Ω")
        print(f"  Voltage at Rectifier Input (Peak) |V_rect_in_{n}| = {V_rect_in_amp:.4f} V")
        print("-" * 25)

    # --- Step 5: Display the final result and the summation equation ---
    equation_str = " + ".join([f"{c:.4f}" for c in contributions])
    print("\n--- Final Calculation ---")
    print("The total voltage is the sum of contributions from each harmonic:")
    print(f"V_total = {equation_str}")
    print(f"V_total = {total_voltage_across_CL:.4f} V")

    return total_voltage_across_CL

# Execute the function and print the final result.
final_voltage = solve_circuit_voltage()
# The final answer in the requested format
# print(f"\n<<< {final_voltage:.4f} >>>") # The calculated value is 1.2729V. Rounding to 2 significant digits as per the source values might be better
print(f"\n<<< {round(final_voltage, 3)} >>>")
