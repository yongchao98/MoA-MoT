import math

def calculate_system_efficiency():
    """
    Calculates the overall system efficiency considering harmonic distortions and parasitic losses.
    """
    # --- System Parameters ---
    V_fundamental = 1.0  # V
    R0_parasitic = 50.0  # Ohms, base parasitic resistance
    f0 = 915e6  # Hz, fundamental frequency
    R_load_eff = 50.0  # Ohms, effective load impedance of the rectifier
    # C_parasitic = 2e-15 # F, This value is not needed for the chosen model, as explained in the plan.
    
    harmonics = [1, 3, 5, 7]
    
    # --- Data Storage ---
    V_amps = []
    R_parasitics = []
    eta_harmonics = []
    P_in_scaled = []
    P_out_scaled = []

    print("--- Calculation Steps ---")
    
    # --- Loop through each harmonic ---
    v_prev = V_fundamental / 0.9 # Initialize for the loop logic
    for n in harmonics:
        # 1. Calculate Harmonic Voltage
        if n == 1:
            v_amp = V_fundamental
        else:
            # Voltage drops by 10% relative to the previous harmonic's amplitude
            v_amp = v_prev * 0.9
        V_amps.append(v_amp)
        v_prev = v_amp

        # 2. Calculate Parasitic Resistance for the harmonic
        # R_p(f) = R0 * (f/f0)^2 = R0 * (n*f0/f0)^2 = R0 * n^2
        r_p = R0_parasitic * (n**2)
        R_parasitics.append(r_p)

        # 3. Calculate Efficiency for the harmonic
        # eta_n = R_load / (R_parasitic + R_load)
        eta_n = R_load_eff / (r_p + R_load_eff)
        eta_harmonics.append(eta_n)
        
        # 4. Calculate Scaled Input and Output Power for weighting
        # P_in is proportional to V^2 / R_total
        p_in_n = (v_amp**2) / (r_p + R_load_eff)
        P_in_scaled.append(p_in_n)
        
        # P_out = P_in * eta
        p_out_n = p_in_n * eta_n
        P_out_scaled.append(p_out_n)

        print(f"\nFor Harmonic n={n}:")
        print(f"  Voltage Amplitude (V_{n}) = {v_amp:.4f} V")
        print(f"  Parasitic Resistance (R_p_{n}) = {r_p:.1f} Ohms")
        print(f"  Harmonic Efficiency (η_{n}) = {eta_n:.4f}")
        print(f"  Scaled Input Power (P_in_{n}) = {p_in_n:.6e}")
        print(f"  Scaled Output Power (P_out_{n}) = {p_out_n:.6e}")

    # 5. Calculate Total Power and Overall Efficiency
    total_P_in = sum(P_in_scaled)
    total_P_out = sum(P_out_scaled)
    overall_efficiency = total_P_out / total_P_in
    
    # --- Final Output ---
    print("\n--- Overall System Efficiency Calculation ---")
    
    p_out_sum_str = " + ".join([f"{p:.6e}" for p in P_out_scaled])
    p_in_sum_str = " + ".join([f"{p:.6e}" for p in P_in_scaled])
    
    print(f"Efficiency = (P_out_1 + P_out_3 + P_out_5 + P_out_7) / (P_in_1 + P_in_3 + P_in_5 + P_in_7)")
    print(f"Efficiency = ({p_out_sum_str}) / ({p_in_sum_str})")
    print(f"Efficiency = {total_P_out:.6e} / {total_P_in:.6e}")
    print(f"\nFinal Overall System Efficiency = {overall_efficiency:.4f} or {overall_efficiency * 100:.2f}%")
    
    # The final answer in the requested format
    print(f"\n<<<{overall_efficiency:.4f}>>>")

calculate_system_efficiency()