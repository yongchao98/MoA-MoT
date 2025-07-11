import numpy as np

def calculate_total_voltage():
    """
    Calculates the total voltage across the load capacitor C_L in the given circuit.

    This function models the system for each harmonic, calculates the power delivered
    by each, sums the power, and then determines the final DC voltage assuming
    100% rectifier efficiency.
    """
    # --- Given Constants ---
    V_RF_peak = 1.0        # V, peak amplitude of the fundamental
    f_fundamental = 915e6  # Hz, fundamental frequency
    R_L = 8000.0           # Ohms
    C_L = 5e-12            # Farads
    
    # --- Parasitic Component Values ---
    R0_parasitic = 50.0    # Ohms, for parasitic resistance calculation
    f0_parasitic = 915e6   # Hz, reference frequency for parasitic resistance
    C_parasitic = 2e-15    # Farads, parasitic capacitance
    
    # --- Derived Constants ---
    # Total capacitance at the load is C_L in parallel with C_parasitic
    C_total_load = C_L + C_parasitic
    
    # Define the harmonics to be considered
    harmonics = [1, 3, 5, 7]
    
    # Calculate the peak voltage of the source for each harmonic
    # Voltage drops by 10% for each higher harmonic (i.e., is 90% of the previous)
    # The exponent is (n-1)/2 because harmonics are odd (1, 3, 5, ...)
    V_source_peaks = [V_RF_peak * (0.9)**((n - 1) // 2) for n in harmonics]

    powers_to_load = []
    
    print("--- Calculating Power Contribution from Each Harmonic ---")
    
    # Loop through each harmonic to calculate its power contribution
    for i, n in enumerate(harmonics):
        # Current harmonic's frequency and peak voltage
        f_n = n * f_fundamental
        V_s_peak = V_source_peaks[i]
        
        # Calculate frequency-dependent parasitic resistance
        R_p = R0_parasitic * (f_n / f0_parasitic)**2
        
        # Calculate the complex impedance of the parallel RC load
        # Z_L = 1 / (1/R_L + j*omega*C_total_load)
        omega_n = 2 * np.pi * f_n
        Z_L_complex = 1 / (1/R_L + 1j * omega_n * C_total_load)
        
        # Calculate the peak voltage across the load using the voltage divider rule
        # V_L = V_s * Z_L / (R_p + Z_L)
        V_L_peak = V_s_peak * abs(Z_L_complex / (R_p + Z_L_complex))
        
        # Convert peak load voltage to RMS voltage
        V_L_rms = V_L_peak / np.sqrt(2)
        
        # Calculate the average power delivered to the load resistor R_L
        # P = V_rms^2 / R_L
        P_n = V_L_rms**2 / R_L
        powers_to_load.append(P_n)
        
        print(f"Harmonic n={n} (f={f_n/1e6:.0f} MHz): Power = {P_n*1e6:.4f} uW")

    # Sum the powers from all harmonics to get the total power delivered to the load
    P_total = sum(powers_to_load)
    
    # Calculate the equivalent DC voltage across the load
    # Assumes 100% AC-to-DC power conversion efficiency: P_total_AC = P_DC = V_DC^2 / R_L
    V_DC = np.sqrt(P_total * R_L)
    
    print("\n--- Final Voltage Calculation ---")
    print(f"Total AC power delivered to R_L: P_total = {P_total*1e6:.4f} uW")
    
    # Display the final calculation in a clear format
    power_sum_str = " + ".join([f"{p:.3e}" for p in powers_to_load])
    print("\nThe total voltage across C_L is calculated as V_DC = sqrt(P_total * R_L)")
    print("Final Equation:")
    print(f"V_total = sqrt( ( {power_sum_str} ) * {R_L:.0f} ) = {V_DC:.4f} V")
    
    return V_DC

if __name__ == '__main__':
    final_voltage = calculate_total_voltage()
    # The final numerical answer is wrapped according to the required format
    print(f"\n<<< {final_voltage:.4f} >>>")
