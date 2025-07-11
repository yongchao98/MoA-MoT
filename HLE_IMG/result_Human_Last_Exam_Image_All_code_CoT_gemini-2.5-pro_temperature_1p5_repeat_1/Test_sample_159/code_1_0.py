import math

def calculate_efficiency():
    """
    Calculates the overall system efficiency considering harmonic distortion and parasitic losses.
    """

    # --- Step 1: Define constants from the problem description ---
    V_1 = 1.0  # V, peak voltage of the fundamental frequency
    harmonic_drop_factor = 0.90  # Voltage drops by 10%, so it's 90% of the previous
    R_0 = 50.0  # Ohms, reference parasitic resistance
    f_0 = 915e6  # Hz, reference frequency for parasitic resistance
    f_1 = 915e6  # Hz, fundamental frequency
    # Assume the rectifier is matched to the system impedance for max power transfer
    R_rect_in = 50.0  # Ohms

    # --- Step 2: Calculate Harmonic Efficiency (eta_harmonic) ---
    print("Step 2: Calculate Harmonic Efficiency\n")
    # Calculate voltage amplitudes of the harmonics
    V_3 = V_1 * harmonic_drop_factor
    V_5 = V_3 * harmonic_drop_factor
    V_7 = V_5 * harmonic_drop_factor

    # Calculate relative power for each harmonic (P is proportional to V^2)
    P_1_rel = V_1**2
    P_3_rel = V_3**2
    P_5_rel = V_5**2
    P_7_rel = V_7**2

    # Calculate total relative power
    P_total_rel = P_1_rel + P_3_rel + P_5_rel + P_7_rel

    # Calculate harmonic efficiency
    eta_harmonic = P_1_rel / P_total_rel
    
    print("Harmonic Efficiency Calculation:")
    print(f"η_harmonic = P1 / (P1 + P3 + P5 + P7)")
    print(f"Relative powers are Pn ∝ Vn^2:")
    print(f"P1 ∝ {V_1:.2f}^2 = {P_1_rel:.4f}")
    print(f"P3 ∝ {V_3:.2f}^2 = {P_3_rel:.4f}")
    print(f"P5 ∝ {V_5:.3f}^2 = {P_5_rel:.4f}")
    print(f"P7 ∝ {V_7:.3f}^2 = {P_7_rel:.4f}")
    print(f"η_harmonic = {P_1_rel:.4f} / ({P_1_rel:.4f} + {P_3_rel:.4f} + {P_5_rel:.4f} + {P_7_rel:.4f})")
    print(f"η_harmonic = {P_1_rel:.4f} / {P_total_rel:.4f} = {eta_harmonic:.4f}\n")


    # --- Step 3: Calculate Parasitic Loss Efficiency (eta_parasitic) ---
    print("Step 3: Calculate Parasitic Loss Efficiency\n")
    # Calculate parasitic resistance at the fundamental frequency f1
    R_parasitic_f1 = R_0 * (f_1 / f_0)**2

    # Calculate parasitic loss efficiency
    eta_parasitic = R_rect_in / (R_parasitic_f1 + R_rect_in)
    
    print("Parasitic Loss Efficiency Calculation:")
    print(f"R_parasitic(f1) = R0 * (f1/f0)^2 = {R_0:.1f} * ({f_1:.0f}/{f_0:.0f})^2 = {R_parasitic_f1:.1f} Ω")
    print("Assuming rectifier input resistance is matched for max power transfer:")
    print(f"R_rect_in = {R_rect_in:.1f} Ω")
    print(f"η_parasitic = R_rect_in / (R_parasitic_f1 + R_rect_in)")
    print(f"η_parasitic = {R_rect_in:.1f} / ({R_parasitic_f1:.1f} + {R_rect_in:.1f}) = {eta_parasitic:.4f}\n")


    # --- Step 4: Calculate Overall System Efficiency ---
    print("Step 4: Calculate Overall System Efficiency\n")
    eta_overall = eta_harmonic * eta_parasitic
    
    print("Overall Efficiency Equation:")
    print("η_overall = η_harmonic * η_parasitic")
    print(f"η_overall = {eta_harmonic:.4f} * {eta_parasitic:.4f}")
    print(f"η_overall = {eta_overall:.4f}")
    print(f"Overall System Efficiency: {eta_overall * 100:.2f}%")
    
    # Final answer in the required format
    print(f"\nFinal calculated efficiency value is {eta_overall}")


# Run the calculation
calculate_efficiency()