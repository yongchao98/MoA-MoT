import math

def calculate_system_efficiency():
    """
    Calculates the overall system efficiency considering harmonic distortions and parasitic losses.
    """
    # --- Parameters from the problem statement ---
    V_rf_peak = 1.0  # V, peak voltage of the fundamental frequency
    f0 = 915e6       # Hz, fundamental frequency
    R0 = 50.0        # Ohms, base parasitic resistance
    # C_parasitic = 2e-15 # Farads, parasitic capacitance (negligible effect)
    
    # Assumption based on problem context
    R_rect = 50.0    # Ohms, assumed input resistance of the rectifier
    
    # Harmonic characteristics
    harmonics = [1, 3, 5, 7]
    voltage_drop_ratio = 0.9 # Voltage drops by 10% for each higher harmonic

    # --- Calculation ---
    
    # Calculate peak voltages for each harmonic
    voltages = {}
    current_v_peak = V_rf_peak
    voltages[1] = current_v_peak
    for i in range(1, len(harmonics)):
        # The voltage of the current harmonic is 0.9 times the voltage of the previous harmonic in the list
        prev_harmonic_v = voltages[harmonics[i-1]]
        current_v_peak = prev_harmonic_v * voltage_drop_ratio
        voltages[harmonics[i]] = current_v_peak

    total_power_input = 0.0
    total_power_delivered_to_rectifier = 0.0
    
    p_rect_list = []
    p_in_list = []

    print("Calculating power components for each harmonic:\n")
    for n in harmonics:
        v_peak = voltages[n]
        
        # Parasitic resistance increases with the square of the harmonic number
        r_parasitic = R0 * (n**2)
        
        # Total resistance in the series model for this harmonic
        r_total = r_parasitic + R_rect
        
        # Peak current for this harmonic
        i_peak = v_peak / r_total
        
        # Power is proportional to I_rms^2 = (I_peak / sqrt(2))^2 = 0.5 * I_peak^2
        # The 0.5 factor will cancel in the final ratio, but we include it for correctness.
        power_factor = 0.5 * (i_peak**2)
        
        # Power delivered to the rectifier at this harmonic
        p_rect_n = power_factor * R_rect
        p_rect_list.append(p_rect_n)
        
        # Total input power from the source at this harmonic
        p_in_n = power_factor * r_total
        p_in_list.append(p_in_n)
        
        total_power_delivered_to_rectifier += p_rect_n
        total_power_input += p_in_n
        
        print(f"Harmonic n={n}:")
        print(f"  V_peak = {v_peak:.3f} V")
        print(f"  R_parasitic = {r_parasitic:.1f} Ω")
        print(f"  R_total = {r_total:.1f} Ω")
        print(f"  P_rectifier(n) = {p_rect_n:.4e} W")
        print(f"  P_input(n) = {p_in_n:.4e} W\n")

    # Calculate overall efficiency
    efficiency = total_power_delivered_to_rectifier / total_power_input

    # --- Final Output ---
    print("--- Overall Efficiency Calculation ---")
    p_rect_str_list = [f"{p:.4e}" for p in p_rect_list]
    print(f"Total power delivered to rectifier P_rect = {' + '.join(p_rect_str_list)}")
    print(f"= {total_power_delivered_to_rectifier:.4e} W")
    
    p_in_str_list = [f"{p:.4e}" for p in p_in_list]
    print(f"\nTotal input power from source P_in = {' + '.join(p_in_str_list)}")
    print(f"= {total_power_input:.4e} W")
    
    print(f"\nOverall Efficiency η = P_rect / P_in")
    print(f"η = {total_power_delivered_to_rectifier:.4e} W / {total_power_input:.4e} W")
    print(f"η = {efficiency:.4f}")
    
    # Return the final answer in the required format
    print(f"\n<<<η = {efficiency:.4f}>>>")


calculate_system_efficiency()