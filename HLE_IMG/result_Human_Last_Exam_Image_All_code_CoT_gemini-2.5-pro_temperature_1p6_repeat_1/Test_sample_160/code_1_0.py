import cmath
import math

def calculate_total_voltage():
    """
    Calculates the total DC voltage across the load capacitor in the given circuit,
    considering harmonic distortion and parasitic effects.
    """

    # --- Given Circuit Parameters ---
    V_RF = 1.0  # V, amplitude of the fundamental harmonic
    f0 = 915e6  # Hz, fundamental frequency
    
    # Parasitic effects parameters
    R0 = 50.0   # Ohm, base parasitic resistance
    C_parasitic = 2e-15  # F (2 fF), parasitic capacitance

    # Circuit components from the table
    Cc = 1e-12  # F (1 pF), coupling capacitor
    RL = 8e3    # Ohm (8 kOhm), load resistor
    # C_L is not needed for the AC analysis under the assumption that RL is the effective load.

    harmonics = [1, 3, 5, 7]
    harmonic_voltages = {}
    
    # --- Step 1: Calculate initial source voltage for each harmonic ---
    v_source = V_RF
    for n in harmonics:
        harmonic_voltages[n] = v_source
        v_source *= 0.9 # Voltage drops by 10% for the next harmonic
    
    print("--- Analysis Setup ---")
    print(f"Fundamental Frequency f0: {f0/1e6} MHz")
    print(f"Load Resistance RL: {RL/1e3} kOhm")
    print(f"Parasitic Resistance R0: {R0} Ohm")
    print(f"Parasitic Capacitance C_parasitic: {C_parasitic*1e15} fF")
    print(f"Coupling Capacitance Cc: {Cc*1e12} pF")
    print("\nSource Voltage Amplitudes:")
    for n, v in harmonic_voltages.items():
        print(f"  - Harmonic {n}: {v:.3f} V")
    print("-" * 24)
    print("\n--- Per-Harmonic Calculation ---")

    total_dc_voltage = 0
    dc_contributions = []

    # --- Step 2: Analyze each harmonic ---
    for n in harmonics:
        Vn = harmonic_voltages[n]
        fn = n * f0
        omega_n = 2 * math.pi * fn

        # Calculate frequency-dependent impedances
        R_p_n = R0 * (n**2)
        
        # Reactance of parasitic capacitance
        Z_Cp_n = 1 / (1j * omega_n * C_parasitic)
        
        # Reactance of two series coupling capacitors
        Z_2Cc_n = 2 / (1j * omega_n * Cc)

        # Model the AC path as a voltage divider.
        # Assumed model: V_n -> R_p_n -> (C_p_n || (Z_2Cc_n + R_L))
        # We want the voltage across R_L.

        # Impedance of the rectifier stage (Z_2Cc + R_L)
        Z_rect_stage = Z_2Cc_n + RL

        # Total load impedance seen by R_p_n is Z_rect_stage in parallel with Z_Cp_n
        Z_total_load = (Z_rect_stage * Z_Cp_n) / (Z_rect_stage + Z_Cp_n)

        # Voltage after the series parasitic resistor R_p_n
        V_after_Rp = Vn * Z_total_load / (R_p_n + Z_total_load)

        # Final voltage across the load resistor R_L
        V_rect_n = V_after_Rp * RL / Z_rect_stage
        
        # The DC contribution is the peak AC amplitude (ideal rectifier)
        V_dc_n = abs(V_rect_n)
        dc_contributions.append(V_dc_n)
        total_dc_voltage += V_dc_n
        
        print(f"For harmonic n={n} (f={fn/1e6:.0f} MHz):")
        print(f"  - Parasitic Resistance R_p = {R_p_n:.1f} Ohm")
        print(f"  - Delivered Peak Voltage to Rectifier = {V_dc_n:.4f} V")

    # --- Step 3: Sum contributions and Print Final Result ---
    print("\n--- Final Calculation ---")
    print("The total voltage across C_L (V_DC) is the sum of contributions from each harmonic.")
    
    equation_str = "V_DC_total = " + " + ".join([f"{v:.4f}" for v in dc_contributions])
    print(f"{equation_str} V")
    print(f"V_DC_total = {total_dc_voltage:.4f} V")

    return total_dc_voltage

if __name__ == "__main__":
    final_voltage = calculate_total_voltage()
    # The final answer format is specified for the platform.
    # To conform, we'll print the numerical answer within the special tags.
    print(f"\n<<<{final_voltage:.4f}>>>")
