import math

def calculate_efficiency():
    """
    Calculates the overall system efficiency considering harmonic distortions and parasitic losses.
    """
    # --- Step 1: Define System Parameters ---
    V_rf = 1.0  # V, peak voltage of the fundamental
    f0 = 915e6  # Hz, fundamental frequency
    R0 = 50.0  # Ohms, base parasitic resistance
    C_parasitic = 2e-15  # F, parasitic capacitance (2 fF)
    R_L = 8e3  # Ohms, load resistance (8 kOhm)
    
    # Assumed rectifier input resistance
    R_rectifier = 50.0 # Ohms

    # Harmonic properties
    harmonics = [1, 3, 5, 7]
    voltages = [
        V_rf,
        V_rf * 0.9,
        V_rf * 0.9**2,
        V_rf * 0.9**3
    ]

    print("--- System Parameters and Assumptions ---")
    print(f"Fundamental Frequency (f0): {f0/1e6} MHz")
    print(f"Load Resistance (R_L): {R_L/1e3} kOhms")
    print(f"Parasitic Resistance (R0): {R0} Ohms")
    print(f"Parasitic Capacitance (C_parasitic): {C_parasitic*1e15} fF")
    print(f"Assumed Rectifier Input Resistance (R_rectifier): {R_rectifier} Ohms")
    print("\nHarmonic Voltages (Peak):")
    for n, v in zip(harmonics, voltages):
        print(f"  V_{n}: {v:.4f} V")
    print("-" * 40)

    # --- Step 2 & 3: Calculate Input and Output Power ---
    total_input_power = 0
    output_power = 0
    
    # Store intermediate results for final equation
    p_in_values = []

    print("\n--- Power Calculation per Harmonic ---")
    for n, v_peak in zip(harmonics, voltages):
        f_n = n * f0
        omega_n = 2 * math.pi * f_n
        
        # Calculate frequency-dependent parasitic resistance
        r_parasitic_n = R0 * (n**2)
        
        # Calculate parasitic capacitive reactance
        x_c_n = -1 / (omega_n * C_parasitic)
        
        # Total series resistance and impedance
        r_total_n = r_parasitic_n + R_rectifier
        z_total_mag_sq = r_total_n**2 + x_c_n**2
        z_total_mag = math.sqrt(z_total_mag_sq)
        
        # RMS voltage and current
        v_rms = v_peak / math.sqrt(2)
        i_rms = v_rms / z_total_mag
        
        # Input power for this harmonic
        p_in_n = i_rms**2 * r_total_n
        total_input_power += p_in_n
        p_in_values.append(p_in_n)
        
        print(f"Harmonic n={n}:")
        print(f"  Frequency: {f_n/1e6:.2f} MHz")
        print(f"  Parasitic Resistance (R_p): {r_parasitic_n:.2f} Ohms")
        print(f"  Parasitic Reactance (X_C): {x_c_n/1e3:.2f} kOhms")
        print(f"  Total Impedance |Z|: {z_total_mag/1e3:.2f} kOhms")
        print(f"  Input Power (P_in_{n}): {p_in_n * 1e9:.4f} nW")

        # Calculate output power (only from fundamental n=1)
        if n == 1:
            v_peak_rect = (v_peak / z_total_mag) * R_rectifier
            v_dc = 2 * v_peak_rect  # Ideal voltage doubler assumption
            output_power = v_dc**2 / R_L
            
            # Store values for final equation
            v_dc_val = v_dc
            v_peak_rect_val = v_peak_rect

    print("-" * 40)
    print("\n--- Total Power and Efficiency Calculation ---")
    
    # --- Step 4: Calculate and Print Final Efficiency ---
    efficiency = output_power / total_input_power

    # Print detailed breakdown of the final calculation
    print(f"Total Input Power (P_in) = P_in_1 + P_in_3 + P_in_5 + P_in_7")
    p_in_str = " + ".join([f"{p*1e6:.4f}" for p in p_in_values])
    print(f"P_in = {p_in_str} uW = {total_input_power * 1e6:.4f} uW\n")

    print(f"Peak Voltage at Rectifier (V_peak_rect,1) = (V_1 / |Z_1|) * R_rectifier")
    print(f"V_peak_rect,1 = ({voltages[0]:.2f} V / {math.sqrt((R0 + R_rectifier)**2 + (-1 / (2 * math.pi * f0 * C_parasitic))**2):.2f} Ohm) * {R_rectifier:.2f} Ohm = {v_peak_rect_val * 1e3:.4f} mV\n")

    print(f"Output DC Voltage (V_DC) = 2 * V_peak_rect,1")
    print(f"V_DC = 2 * {v_peak_rect_val * 1e3:.4f} mV = {v_dc_val * 1e3:.4f} mV\n")

    print(f"Output Power (P_out) = V_DC^2 / R_L")
    print(f"P_out = ({v_dc_val:.6f} V)^2 / {R_L:.0f} Ohm = {output_power * 1e9:.4f} nW\n")

    print("Overall Efficiency (η) = P_out / P_in")
    print(f"η = {output_power:.4e} W / {total_input_power:.4e} W")
    print(f"η = {efficiency:.6f}")
    print(f"η = {efficiency * 100:.6f} %")
    
    return efficiency * 100

# Execute the calculation
final_efficiency_percent = calculate_efficiency()
# The final answer is requested in a specific format.
# print(f"\n<<<>>>")
# The final answer is the efficiency in percent.
# For example, if the efficiency is 0.002863%, the answer is 0.002863
final_answer = f"{final_efficiency_percent:.4f}"
# print(f"\n<<<{final_answer}>>>")