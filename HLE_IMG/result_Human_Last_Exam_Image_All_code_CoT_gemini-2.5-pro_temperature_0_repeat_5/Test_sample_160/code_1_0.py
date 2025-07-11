import numpy as np

def calculate_total_voltage():
    """
    Calculates the total RMS voltage across the load capacitor considering
    harmonics and parasitic effects.
    """
    # --- Parameters ---
    # From the problem description and image table
    V_rf = 1.0          # V, peak voltage of the fundamental
    f1 = 915e6          # Hz, fundamental frequency
    R_L = 8e3           # Ohms
    C_L = 5e-12         # Farads
    R0 = 50.0           # Ohms, parasitic resistance base value
    f0 = 915e6          # Hz, parasitic resistance reference frequency
    C_parasitic = 2e-15 # Farads, parasitic capacitance

    # --- Calculations ---

    # Total capacitance in parallel with the load resistor
    C_total = C_L + C_parasitic

    harmonics = [1, 3, 5, 7]
    V_L_peak_sq_list = []

    print("--- Calculation Breakdown ---")
    print(f"Load Resistor R_L = {R_L:.0f} Ohms")
    print(f"Load Capacitor C_L = {C_L*1e12:.1f} pF")
    print(f"Parasitic Capacitor C_parasitic = {C_parasitic*1e15:.1f} fF")
    print(f"Total Load Capacitance C_total = {C_total*1e12:.4f} pF\n")

    # Loop through each harmonic
    for n in harmonics:
        # Frequency of the current harmonic
        f_n = n * f1
        omega_n = 2 * np.pi * f_n

        # Source voltage amplitude for the current harmonic
        # V_n = V_1 * (0.9)^((n-1)/2)
        V_S_n = V_rf * (0.9)**((n - 1) / 2)

        # Frequency-dependent parasitic resistance
        R_p_n = R0 * (f_n / f0)**2

        # Load impedance for the current harmonic: Z_L = R_L || (1 / jw_n C_total)
        Z_L_n = 1 / (1/R_L + 1j * omega_n * C_total)

        # Voltage across the load using voltage divider rule
        V_L_n = V_S_n * Z_L_n / (Z_L_n + R_p_n)

        # Get the peak voltage (magnitude)
        V_L_n_peak = np.abs(V_L_n)

        # Store the square of the peak voltage
        V_L_peak_sq_list.append(V_L_n_peak**2)

        print(f"--- Harmonic n = {n} ---")
        print(f"Frequency f_{n} = {f_n/1e6:.0f} MHz")
        print(f"Source Voltage V_S,{n} (peak) = {V_S_n:.4f} V")
        print(f"Parasitic Resistance R_p,{n} = {R_p_n:.2f} Ohms")
        print(f"Load Impedance |Z_L,{n}| = {np.abs(Z_L_n):.4f} Ohms")
        print(f"Load Voltage V_L,{n} (peak) = {V_L_n_peak:.4f} V")
        print(f"(V_L,{n}_peak)^2 = {V_L_n_peak**2:.6f} V^2\n")

    # Sum of the squared peak voltages
    V_L_peak_sq_sum = sum(V_L_peak_sq_list)

    # The total voltage is the RMS value of the sum of the harmonic components
    V_total_rms = np.sqrt(0.5 * V_L_peak_sq_sum)

    print("--- Final Calculation ---")
    print("The total RMS voltage is the root of the sum of the mean squares of each harmonic component:")
    print("V_total_rms = sqrt( 0.5 * [ (V_L,1_peak)^2 + (V_L,3_peak)^2 + (V_L,5_peak)^2 + (V_L,7_peak)^2 ] )")
    print(f"V_total_rms = sqrt( 0.5 * [ {V_L_peak_sq_list[0]:.6f} + {V_L_peak_sq_list[1]:.6f} + {V_L_peak_sq_list[2]:.6f} + {V_L_peak_sq_list[3]:.6f} ] )")
    print(f"V_total_rms = sqrt( 0.5 * {V_L_peak_sq_sum:.6f} )")
    print(f"V_total_rms = sqrt( {0.5 * V_L_peak_sq_sum:.6f} )")
    print(f"\nFinal Answer: The total voltage across C_L is {V_total_rms:.4f} V (RMS).")
    print(f"<<<{V_total_rms:.4f}>>>")

if __name__ == '__main__':
    calculate_total_voltage()