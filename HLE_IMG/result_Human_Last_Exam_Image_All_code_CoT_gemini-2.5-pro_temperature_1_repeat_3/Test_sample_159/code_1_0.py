import math

def calculate_efficiency():
    """
    Calculates the overall system efficiency considering harmonic distortions and parasitic losses.
    """
    # --- Problem Parameters from the user and the image ---
    V_rf_fund = 1.0     # V, amplitude of the fundamental
    f_fund = 915e6      # Hz, fundamental frequency
    R0 = 50.0           # Ohm, base parasitic resistance
    f0 = 915e6          # Hz, reference frequency for R_parasitic
    C_parasitic = 2e-15 # F (2 fF), parasitic capacitance
    R_L = 8000.0        # Ohm (8 kOhm), load resistance

    harmonics = [1, 3, 5, 7]
    voltage_drop_factor = 0.9 # Voltage drops by 10% for each higher harmonic

    # --- Engineering Assumption ---
    # The AC input impedance of the rectifier stage is not provided.
    # A common approximation for a full-wave rectifier is that its effective
    # input resistance is half the DC load resistance.
    R_in_rectifier = R_L / 2
    print(f"--- ASSUMPTION ---")
    print(f"Assuming the rectifier's effective AC input resistance R_in_rectifier = R_L / 2 = {R_in_rectifier:.0f} Ω\n")

    # --- Calculation ---
    total_input_power = 0
    total_rectifier_power = 0

    print("--- Calculating Power for Each Harmonic ---")
    for n in harmonics:
        # 1. Calculate harmonic frequency and voltage amplitude
        f_n = n * f_fund
        # The voltage drops by 10% for each *higher* harmonic, meaning V3=0.9*V1, V5=0.9*V3, etc.
        # This corresponds to a power series for the harmonic steps (1->3 is 1 step, 1->5 is 2 steps).
        harmonic_step = (n - 1) / 2
        V_n = V_rf_fund * (voltage_drop_factor ** harmonic_step)

        # 2. Calculate frequency-dependent parasitics
        R_p = R0 * (f_n / f0)**2
        # Reactance of parasitic capacitor: Xc = 1 / (omega * C) = 1 / (2*pi*f*C)
        Xc_p = 1 / (2 * math.pi * f_n * C_parasitic)

        # 3. Calculate total series impedance for this harmonic
        # Z_in = R_in_rectifier + R_p + Z_c = (R_in_rectifier + R_p) - j*Xc_p
        Z_real = R_p + R_in_rectifier
        Z_imag = -Xc_p
        Z_mag_sq = Z_real**2 + Z_imag**2

        # 4. Calculate power based on peak voltage and impedance
        # P = V_rms^2 / R_eff = (V_peak / sqrt(2))^2 / R_eff = V_peak^2 / (2 * R_eff)
        # For a complex impedance Z, the current is I = V / Z.
        # Power dissipated in a resistive component R is P = I_rms^2 * R = |I_peak|^2/2 * R
        # |I_peak|^2 = |V_peak / Z|^2 = V_peak^2 / |Z|^2
        # So P = (V_n**2 / Z_mag_sq / 2) * R
        
        P_in_n = (V_n**2 / (2 * Z_mag_sq)) * Z_real           # Total power drawn from source
        P_rect_n = (V_n**2 / (2 * Z_mag_sq)) * R_in_rectifier # Power delivered to rectifier

        # 5. Accumulate total powers
        total_input_power += P_in_n
        total_rectifier_power += P_rect_n
        
        print(f"Harmonic n={n}: V={V_n:.3f}V, f={f_n/1e6:.1f}MHz")
        print(f"  P_input = {P_in_n * 1e6:.4f} µW")
        print(f"  P_delivered_to_rectifier = {P_rect_n * 1e6:.4f} µW\n")


    # --- Final Efficiency Calculation ---
    # The overall system efficiency is the ratio of the power that reaches the
    # ideal rectifier to the total power drawn from the source.
    system_efficiency = total_rectifier_power / total_input_power

    print("--- Total Power Summary ---")
    print(f"Total Input Power (from source): {total_input_power * 1e6:.4f} µW")
    print(f"Total Power Delivered to Rectifier: {total_rectifier_power * 1e6:.4f} µW\n")
    
    print("--- Final Result ---")
    print(f"Overall System Efficiency = Total Power Delivered to Rectifier / Total Input Power")
    print(f"Efficiency = {total_rectifier_power:.4e} W / {total_input_power:.4e} W")
    print(f"Efficiency = {system_efficiency:.4f}")

    return system_efficiency

if __name__ == '__main__':
    efficiency = calculate_efficiency()
    print(f"\n<<<The calculated efficiency is approximately {efficiency*100:.2f}%.>>>")
    print(f"<<<{efficiency:.4f}>>>")
