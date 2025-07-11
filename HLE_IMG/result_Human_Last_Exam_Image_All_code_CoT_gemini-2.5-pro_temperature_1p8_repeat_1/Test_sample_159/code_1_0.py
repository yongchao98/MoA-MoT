import numpy as np

def calculate_system_efficiency():
    """
    Calculates the overall system efficiency considering harmonic distortions and parasitic losses.
    """
    # --- 1. Define Constants and Parameters ---
    f0 = 915e6  # Fundamental frequency in Hz
    V_rf = 1.0  # Fundamental peak voltage in Volts
    
    # Parasitic and harmonic distortion parameters from the text
    R0 = 50.0  # Parasitic resistance base value in Ohms
    f_R0 = 915e6 # Reference frequency for R_parasitic
    C_parasitic = 2e-15  # Parasitic capacitance in Farads
    voltage_drop_factor = 0.9  # Voltage of next harmonic is 90% of the previous one

    # Assumption: The rectifier input impedance is a 50 Ohm resistor.
    # This is based on the hint that R_parasitic(f0) = R0 = 50 Ohms, suggesting a 50 Ohm system.
    R_rectifier_in = 50.0

    # --- 2. Setup Harmonics ---
    harmonics = [1, 3, 5, 7]
    voltages = {}
    voltages[1] = V_rf
    for i in range(1, len(harmonics)):
        prev_n = harmonics[i-1]
        current_n = harmonics[i]
        voltages[current_n] = voltages[prev_n] * voltage_drop_factor

    total_proportional_power_out = 0
    total_proportional_power_in = 0

    print("--- System Parameters ---")
    print(f"Fundamental Frequency (f0): {f0 / 1e6} MHz")
    print(f"Parasitic Resistance (R0 @ f0): {R0} Ohms")
    print(f"Parasitic Capacitance (C_parasitic): {C_parasitic * 1e15} fF")
    print(f"Assumed Rectifier Input Resistance (R_rectifier_in): {R_rectifier_in} Ohms")
    print("-" * 40)

    # --- 3. Loop Through Harmonics and Calculate Powers ---
    for n in harmonics:
        print(f"--- Calculations for Harmonic n={n} ---")
        f_n = n * f0
        omega_n = 2 * np.pi * f_n
        V_n = voltages[n]

        # Calculate frequency-dependent parasitic resistance
        R_p_n = R0 * (f_n / f_R0)**2

        # Calculate the effective input resistance of the rectifier, considering the shunt parasitic capacitance.
        # This is the real part of the parallel combination of R_rectifier_in and C_parasitic.
        omega_C_R = omega_n * C_parasitic * R_rectifier_in
        R_L_eff_n = R_rectifier_in / (1 + omega_C_R**2)

        # Calculate the total impedance Z_total = R_p_n + (R_rectifier_in || C_parasitic)
        Z_total_real = R_p_n + R_L_eff_n
        Z_total_imag = -R_L_eff_n * omega_C_R
        Z_total_mod_sq = Z_total_real**2 + Z_total_imag**2

        # Calculate proportional power delivered to the rectifier's resistance.
        # This is our "output" power for the network. P_out is proportional to V_n^2 * R_L_eff_n / |Z_total_n|^2
        power_out_n_term = V_n**2 * R_L_eff_n / Z_total_mod_sq
        total_proportional_power_out += power_out_n_term

        # Calculate proportional power drawn from the source.
        # This is our "input" power for the network. P_in is proportional to V_n^2 * Re{Z_total_n} / |Z_total_n|^2
        power_in_n_term = V_n**2 * Z_total_real / Z_total_mod_sq
        total_proportional_power_in += power_in_n_term
        
        print(f"  Frequency (f_n): {f_n / 1e6:.2f} MHz, Voltage (V_n): {V_n:.3f} V")
        print(f"  Parasitic R (R_p_n): {R_p_n:.2f} Ohms")
        print(f"  Effective Load R (R_L_eff_n): {R_L_eff_n:.2f} Ohms")
        print(f"  Total Z Real Part: {Z_total_real:.2f} Ohms")
        print("-" * 40)


    # --- 4. Calculate Final Efficiency ---
    if total_proportional_power_in == 0:
        efficiency = 0
    else:
        efficiency = total_proportional_power_out / total_proportional_power_in

    print("\n--- Overall Efficiency Calculation ---")
    print(f"Total Proportional Power Delivered to Rectifier (Sum of P_out terms): {total_proportional_power_out:.6f}")
    print(f"Total Proportional Power Drawn from Source (Sum of P_in terms): {total_proportional_power_in:.6f}")
    
    print("\nFinal Equation:")
    print(f"Efficiency = Total Delivered Power / Total Input Power")
    print(f"Efficiency = {total_proportional_power_out:.6f} / {total_proportional_power_in:.6f}")
    
    print(f"\nResult:")
    print(f"Calculated Overall System Efficiency = {efficiency:.4f}")
    print(f"Efficiency in Percent = {efficiency * 100:.1f}%")

    return efficiency * 100

# Execute the calculation
final_efficiency_percent = calculate_system_efficiency()
# print(f"<<<{final_efficiency_percent:.1f}>>>")