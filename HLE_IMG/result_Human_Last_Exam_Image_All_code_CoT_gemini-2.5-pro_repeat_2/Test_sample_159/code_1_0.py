import math

def calculate_efficiency():
    """
    Calculates the overall system efficiency considering harmonic distortions and parasitic losses.
    """
    # --- Step 1: Define System Constants ---
    V_rf = 1.0  # Peak voltage of the fundamental (V)
    f0 = 915e6  # Fundamental frequency (Hz)
    R_L = 8000.0  # Load resistance, used to approximate rectifier input resistance (Ohm)
    R0_parasitic = 50.0  # Reference parasitic resistance (Ohm)
    f0_parasitic = 915e6 # Reference frequency for R_parasitic (Hz)
    C_parasitic = 2e-15  # Parasitic capacitance (F)
    
    w0 = 2 * math.pi * f0

    # --- Step 2: Model Harmonic Voltages ---
    harmonics = [1, 3, 5, 7]
    voltages = {1: V_rf}
    # The voltage drops by 10% for each higher harmonic
    for i in range(1, len(harmonics)):
        prev_harmonic = harmonics[i-1]
        current_harmonic = harmonics[i]
        voltages[current_harmonic] = voltages[prev_harmonic] * 0.9

    total_input_power = 0
    delivered_power_at_fundamental = 0

    print("Calculating efficiency based on power transfer to the rectifier input.\n")

    # --- Step 3 & 4: Calculate Impedance and Power for Each Harmonic ---
    for n in harmonics:
        V_n = voltages[n]
        f_n = n * f0
        w_n = n * w0

        # Parasitic resistance at harmonic frequency
        R_p_n = R0_parasitic * (f_n / f0_parasitic)**2

        # Impedance of the load model (R_L || C_parasitic)
        # Admittance Y_load = 1/R_L + j*w*C
        # Impedance Z_load = 1 / Y_load
        # A more direct calculation for the real and imaginary parts of Z_load:
        # Z_load = R_L / (1 + j*w*C*R_L)
        wcR = w_n * C_parasitic * R_L
        Re_Z_load_n = R_L / (1 + wcR**2)
        Im_Z_load_n = -wcR * R_L / (1 + wcR**2)

        # Total input impedance (R_parasitic in series with Z_load)
        Re_Z_total_n = R_p_n + Re_Z_load_n
        Im_Z_total_n = Im_Z_load_n
        mag_sq_Z_total_n = Re_Z_total_n**2 + Im_Z_total_n**2

        # Total power from the source at this harmonic
        # P_in = (V_peak^2 / 2) * Re(Z_total) / |Z_total|^2
        P_in_n = (V_n**2 / 2) * Re_Z_total_n / mag_sq_Z_total_n
        total_input_power += P_in_n
        
        # --- Step 5: Calculate "Useful" Power (at n=1) ---
        if n == 1:
            # Power delivered to the rectifier load model (Z_load) at the fundamental frequency
            # P_delivered = |I|^2 * Re(Z_load) = (V_peak^2 / (2*|Z_total|^2)) * Re(Z_load)
            delivered_power_at_fundamental = (V_n**2 / 2) * Re_Z_load_n / mag_sq_Z_total_n

    # --- Step 6: Calculate Overall Efficiency ---
    if total_input_power == 0:
        efficiency = 0
    else:
        efficiency = delivered_power_at_fundamental / total_input_power

    print("--- Power Summary ---")
    print(f"Total power drawn from source (all harmonics): {total_input_power * 1e6:.2f} uW")
    print(f"Power delivered to rectifier at fundamental: {delivered_power_at_fundamental * 1e6:.2f} uW\n")
    
    print("--- Final Efficiency Calculation ---")
    print("Efficiency = Power delivered at fundamental / Total input power")
    print(f"Efficiency = {delivered_power_at_fundamental:.4g} W / {total_input_power:.4g} W")
    print(f"The overall system efficiency is {efficiency:.4f}, or {efficiency*100:.2f}%.")
    
    return efficiency

# Execute the calculation
final_efficiency = calculate_efficiency()
print(f"\n<<<__{final_efficiency:.4f}__>>>")