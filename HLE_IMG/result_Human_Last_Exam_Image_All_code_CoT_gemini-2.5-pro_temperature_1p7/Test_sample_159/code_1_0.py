import math

def calculate_efficiency():
    # --- Step 1: Define Constants and Assumptions ---
    V_RF = 1.0  # V, amplitude of the fundamental
    f0 = 915e6  # Hz, fundamental frequency
    R0 = 50.0   # Ohm, base parasitic resistance
    # Assume the rectifier input impedance is matched to R0
    R_in_rect = 50.0 # Ohm

    # --- Step 2: Calculate Harmonic Properties ---
    harmonics = [1, 3, 5, 7]
    voltages = {}
    parasitic_resistances = {}
    total_impedances = {}
    
    # Calculate voltage amplitudes based on the 10% drop rule
    v_current = V_RF
    voltages[1] = v_current
    for i in range(len(harmonics) - 1):
        v_next = v_current * 0.9
        voltages[harmonics[i+1]] = v_next
        v_current = v_next

    # Calculate resistances and impedances for each harmonic
    for n in harmonics:
        f_n = n * f0
        # R_parasitic(f) = R0 * (f / f0)^2
        R_p_n = R0 * (f_n / f0)**2
        parasitic_resistances[n] = R_p_n
        total_impedances[n] = R_p_n + R_in_rect

    # --- Step 3: Calculate Power Components ---
    # The actual power is P = 0.5 * V^2 / Z for real impedance.
    # The constant 0.5 cancels out in the efficiency ratio, so we can work with relative powers P' = V^2 / Z.
    
    source_powers = {}
    for n in harmonics:
        P_source_n = voltages[n]**2 / total_impedances[n]
        source_powers[n] = P_source_n
        
    total_source_power = sum(source_powers.values())

    # Calculate useful power delivered to the rectifier at the fundamental frequency (n=1)
    # P_rect_in = 0.5 * I^2 * R_in = 0.5 * (V/Z_total)^2 * R_in
    # Relative power P'_rect_in = (V^2 / Z_total^2) * R_in
    n1 = 1
    rectifier_power_in_f1 = (voltages[n1]**2 / total_impedances[n1]**2) * R_in_rect

    # --- Step 4: Compute Final Efficiency and Print Results ---
    efficiency = rectifier_power_in_f1 / total_source_power

    print("--- Calculation of Overall System Efficiency ---")
    print("\nAssumptions:")
    print(f"  - Rectifier Input Impedance (R_in): {R_in_rect} Ohm")
    print("  - Parasitic Capacitance is negligible.")
    print("  - Rectifier Power Conversion Efficiency is ideal (100%).")

    print("\nHarmonic Analysis:")
    for n in harmonics:
        print(f"  Harmonic n={n}:")
        print(f"    - Voltage (V{n}): {voltages[n]:.3f} V")
        print(f"    - Parasitic Resistance (R_p{n}): {parasitic_resistances[n]:.1f} Ohm")
        print(f"    - Total Impedance (Z{n} = R_p{n} + R_in): {total_impedances[n]:.1f} Ohm")

    print("\nPower Calculation (relative units):")
    print("Total Source Power = P_source(f1) + P_source(f3) + P_source(f5) + P_source(f7)")
    p_eq = " + ".join([f"{p:.6f}" for p in source_powers.values()])
    print(f"Total Source Power = {p_eq} = {total_source_power:.6f}")

    print("\nUseful Power Delivered to Rectifier (at fundamental):")
    print(f"P_rect_in(f1) = (V1^2 / Z1^2) * R_in = ({voltages[1]:.1f}^2 / {total_impedances[1]:.1f}^2) * {R_in_rect:.1f} = {rectifier_power_in_f1:.6f}")
    
    print("\nOverall System Efficiency (η):")
    print(f"η = P_rect_in(f1) / Total Source Power")
    print(f"η = {rectifier_power_in_f1:.6f} / {total_source_power:.6f}")
    print(f"η = {efficiency:.4f} or {efficiency:.2%}")

    # Return the final numerical answer as a percentage
    return efficiency * 100

final_efficiency = calculate_efficiency()
# <<<40.53>>>