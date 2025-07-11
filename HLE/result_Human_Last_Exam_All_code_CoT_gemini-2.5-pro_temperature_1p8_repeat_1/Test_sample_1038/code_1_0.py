import math

# This script demonstrates the quantitative advantage of using on-chip body biasing (Strategy C)
# to improve the operating margin of a transistor in a low-voltage circuit.

def analyze_margins():
    """
    Calculates and compares the safety margin for a transistor to stay in the subthreshold
    region, with and without body biasing, when subjected to a voltage offset.
    """
    # --- Givens and Assumptions based on the problem description ---
    # All voltages are in Volts.
    vt0 = 0.45          # Original Threshold Voltage of the transistor.
    v_offset = 0.1      # Sensor offset is +/- 100 mV.

    # Assumed typical body effect parameters for a standard CMOS process.
    gamma = 0.4         # Body effect coefficient (V^0.5).
    phi_f_2 = 0.7       # 2 * Fermi Potential (V).

    # From Strategy C, a reverse body bias is applied.
    # For an NMOS, this is a positive Source-to-Body voltage (V_sb).
    v_sb_bias = 0.3     # Applied source-to-body reverse bias.

    # We assume a central Vgs operating point is chosen. Due to the +/-100mV offset,
    # the effective Vgs can swing. To maintain high resistance, Vgs must *always*
    # be less than Vt. We'll analyze the worst-case Vgs (Vgs_max).
    # Let's assume a design where the maximum Vgs under offset conditions reaches 0.4V.
    vgs_max = 0.40

    # --- Scenario 1: Standard Transistor (No Body Bias) ---
    print("--- Analysis of Standard Transistor (without Strategy C) ---")
    vt_no_bias = vt0
    # The safety margin is the difference between the threshold and the peak Vgs.
    # A small or negative margin indicates high risk of leaving subthreshold operation.
    margin_no_bias = vt_no_bias - vgs_max
    
    print(f"Original Threshold Voltage (Vt): {vt_no_bias:.3f} V")
    print(f"Maximum Gate-Source Voltage (Vgs_max) due to offset: {vgs_max:.3f} V")
    print(f"Safety Margin (Vt - Vgs_max): {vt_no_bias:.3f} V - {vgs_max:.3f} V = {margin_no_bias:.3f} V")
    if margin_no_bias <= 0.05:
        print("Result: The margin is dangerously small. The transistor is very close to or in strong inversion, compromising high resistance.")
    else:
        print("Result: The transistor remains in subthreshold, but the margin is tight.")
    print("-" * 60)

    # --- Scenario 2: With On-Chip Body Bias (Strategy C) ---
    print("--- Analysis of Transistor with Body Biasing (Strategy C) ---")
    # Equation for new threshold voltage: Vt_new = Vt0 + gamma * (sqrt(2*phi_f + V_sb) - sqrt(2*phi_f))
    sqrt_term_1 = math.sqrt(phi_f_2 + v_sb_bias)
    sqrt_term_2 = math.sqrt(phi_f_2)
    vt_new = vt0 + gamma * (sqrt_term_1 - sqrt_term_2)

    # The new safety margin with the elevated threshold voltage.
    margin_with_bias = vt_new - vgs_max

    print("The equation for the new threshold voltage is: Vt_new = Vt0 + gamma * (sqrt(2*phi_f + V_sb) - sqrt(2*phi_f))")
    print(f"Calculating Vt_new: {vt0:.3f} + {gamma:.2f} * (sqrt({phi_f_2:.2f} + {v_sb_bias:.2f}) - sqrt({phi_f_2:.2f})) = {vt_new:.3f} V")
    print(f"New Elevated Threshold Voltage (Vt_new): {vt_new:.3f} V")
    print(f"Maximum Gate-Source Voltage (Vgs_max) remains: {vgs_max:.3f} V")
    print(f"New Safety Margin (Vt_new - Vgs_max): {vt_new:.3f} V - {vgs_max:.3f} V = {margin_with_bias:.3f} V")
    if margin_with_bias > 0.1:
         print("Result: The margin is now significantly larger, ensuring robust subthreshold operation.")
    else:
         print("Result: Margin has improved.")

# Run the analysis
analyze_margins()