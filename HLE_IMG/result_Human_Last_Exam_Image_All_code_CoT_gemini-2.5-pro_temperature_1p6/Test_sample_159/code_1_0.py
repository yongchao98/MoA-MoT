import math

def calculate_efficiency():
    """
    Calculates the overall system efficiency considering harmonic distortions and parasitic losses.
    """

    # --- Initial Parameters ---
    V_RF = 1.0  # V, amplitude of the fundamental frequency
    f0 = 915e6  # Hz, fundamental frequency
    R0 = 50.0   # Ohms, base parasitic resistance
    # Assume the rectifier input is matched and presents a 50 Ohm load.
    R_rect = 50.0 # Ohms
    
    harmonics = [1, 3, 5, 7]
    voltages = {}
    frequencies = {}
    R_parasitic = {}
    
    # --- Calculate Harmonic Voltages ---
    # Voltage drops by 10% for each higher harmonic
    v_prev = V_RF / 0.9 # Initialize for the loop
    for n in harmonics:
        v_current = v_prev * 0.9
        voltages[n] = v_current
        v_prev = v_current

    # --- Calculate Parameters for each Harmonic ---
    for n in harmonics:
        # Frequency
        frequencies[n] = n * f0
        # Parasitic Resistance R_parasitic(f) = R0 * (f/f0)^2
        R_parasitic[n] = R0 * (frequencies[n] / f0)**2

    # --- Calculate Power for each Harmonic ---
    P_in = {}
    P_rect = {}
    
    for n in harmonics:
        # Total series resistance Z_n = R_parasitic_n + R_rect
        Z_n = R_parasitic[n] + R_rect
        # Total input power P_in = V_rms^2 / Z_n = (V_peak^2 / 2) / Z_n
        P_in[n] = (voltages[n]**2 / 2) / Z_n
        # Power delivered to rectifier (using power divider rule)
        # P_rect = P_in * R_rect / (R_parasitic + R_rect)
        P_rect[n] = P_in[n] * (R_rect / Z_n)
        
    # --- Calculate Total Power and Efficiency ---
    P_in_total = sum(P_in.values())
    P_rect_total = sum(P_rect.values())
    
    efficiency = P_rect_total / P_in_total
    
    # --- Print the Results ---
    print("--- Calculation of System Efficiency ---\n")
    print(f"Assumed rectifier input resistance: {R_rect} Ohms\n")
    
    print("Harmonic | Freq (MHz) | Voltage (V) | R_parasitic (Ohm) | P_input (W)    | P_rectifier (W)")
    print("-----------------------------------------------------------------------------------------")
    for n in harmonics:
        print(f"{n:<8} | {frequencies[n]/1e6:<10.0f} | {voltages[n]:<11.1f} | {R_parasitic[n]:<17.1f} | {P_in[n]:<14.3e} | {P_rect[n]:.3e}")
        
    print("\n--- Final Efficiency Calculation ---\n")
    
    p_rect_sum_str = " + ".join([f"{p:.3e}" for p in P_rect.values()])
    print(f"Total Power Delivered to Rectifier (P_rect_total):")
    print(f"{p_rect_sum_str} = {P_rect_total:.3e} W\n")
    
    p_in_sum_str = " + ".join([f"{p:.3e}" for p in P_in.values()])
    print(f"Total Input Power from Source (P_in_total):")
    print(f"{p_in_sum_str} = {P_in_total:.3e} W\n")
    
    print("Overall System Efficiency (η) = P_rect_total / P_in_total")
    print(f"η = {P_rect_total:.3e} W / {P_in_total:.3e} W")
    print(f"η = {efficiency:.4f}")
    print(f"η = {efficiency * 100:.2f} %")

    # Final answer for the platform
    # print(f"\n<<<{efficiency * 100:.1f}>>>")

calculate_efficiency()