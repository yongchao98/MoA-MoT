import numpy as np

def solve_circuit_voltage():
    """
    Calculates the total DC voltage across the load capacitor in the given circuit,
    considering harmonic distortion and parasitic effects.
    """
    # --- Constants from the problem description and circuit diagram ---
    V_RF = 1.0         # V (Fundamental input voltage amplitude)
    f0 = 915e6         # Hz (Fundamental frequency)
    R0 = 50.0          # Ohms (Base parasitic resistance)
    C_parasitic = 2e-15 # F (2 fF, Parasitic capacitance)
    Cc = 1e-12         # F (1 pF, Coupling capacitance)
    RL = 8000.0        # Ohms (Load resistance)
    CL = 5e-12         # F (5 pF, Load capacitance)
    
    # --- Calculation ---
    harmonics = [1, 3, 5, 7]
    V_dc_total = 0.0
    V_load_harmonics = []

    print("Calculating the total voltage across the load capacitor C_L.")
    print("This is modeled as the sum of the rectified DC voltages from each harmonic component.")
    print("V_DC_total = V_DC(f1) + V_DC(f3) + V_DC(f5) + V_DC(f7)\n")
    
    # The parasitic capacitance is in parallel with the load capacitor
    C_total_load = CL + C_parasitic

    for n in harmonics:
        # --- Per-harmonic calculations ---
        f_n = n * f0
        omega_n = 2 * np.pi * f_n
        
        # Calculate the input voltage for this harmonic, with a 10% drop for each step
        harmonic_step = (n - 1) / 2
        V_n = V_RF * (0.9 ** harmonic_step)
        
        # Calculate the frequency-dependent parasitic resistance
        R_p_n = R0 * (n**2)
        
        # --- Calculate complex impedances ---
        # Load impedance: RL || (CL + C_parasitic)
        Z_C_load_n = 1 / (1j * omega_n * C_total_load)
        Z_load_n = (RL * Z_C_load_n) / (RL + Z_C_load_n)
        
        # Series coupling capacitors impedance (two Cc capacitors in series)
        Z_cc_series_n = 2 / (1j * omega_n * Cc)
        
        # Total impedance of the AC path (voltage divider)
        Z_total_n = R_p_n + Z_cc_series_n + Z_load_n
        
        # Calculate the AC voltage amplitude across the load using the voltage divider rule.
        # This is assumed to be the contribution to the DC voltage after ideal rectification.
        V_load_n = V_n * abs(Z_load_n) / abs(Z_total_n)
        V_load_harmonics.append(V_load_n)
        V_dc_total += V_load_n
        
        # --- Print step-by-step results for this harmonic ---
        print(f"--- Harmonic n={n} ---")
        print(f"Frequency f_{n} = {f_n/1e6:.2f} MHz")
        print(f"Input Voltage V_{n} = {V_n:.4f} V")
        print(f"Parasitic Resistance R_p({n}) = {R_p_n:.2f} Ohms")
        print(f"Load Impedance |Z_load({n})| = {abs(Z_load_n):.4f} Ohms")
        print(f"Total Path Impedance |Z_total({n})| = {abs(Z_total_n):.4f} Ohms")
        print(f"Voltage contribution V_DC({n}) = {V_n:.4f} V * {abs(Z_load_n):.4f} Ohms / {abs(Z_total_n):.4f} Ohms = {V_load_n:.4f} V\n")

    # --- Final Result ---
    equation_str = " + ".join([f"{v:.4f}" for v in V_load_harmonics])
    print("Final Calculation:")
    print(f"V_DC_total = {equation_str}")
    print(f"V_DC_total = {V_dc_total:.4f} V")
    
    # Print the final answer in the required format
    print(f"\n<<<{V_dc_total:.4f}>>>")

solve_circuit_voltage()