import math

def analyze_design_strategies():
    """
    Analyzes and compares different design strategies for a bootstrapped
    pseudo-resistor by modeling key performance metrics.
    """
    # --- Parameters ---
    # These parameters are based on the problem description and standard physics.
    n = 1.5  # Subthreshold slope factor (given range 1.4-1.6)
    V_T = 0.026  # Thermal voltage in Volts (kT/q at room temp)
    n_V_T = n * V_T  # The subthreshold ideality term in Volts
    sensor_offset = 0.100  # 100 mV offset
    charge_injection_q = 10e-15 # Assumed injected charge in Coulombs (10 fC), a typical value
    cap_A = 1e-12 # Capacitor for strategy A in Farads (1 pF)
    # For strategy E, the efficient split-gate reset allows a larger capacitor
    # to be used without violating the reset time constraint.
    cap_E = 10e-12 # Assumed larger capacitor for strategy E (10 pF)

    print("--- Analysis of Pseudo-Resistor Design Strategies ---")
    print(f"Assumed parameters: n={n}, V_T={V_T}V, n*V_T={n_V_T:.3f}V\n")

    # --- Step 1: Analyze Offset Tolerance (Comparing D to Bootstrapped designs) ---
    print("--- Step 1: Impact of 100mV Sensor Offset on Resistance ---")
    print("Strategy D uses a fixed gate bias. A sensor offset on the source node (Vs) directly changes Vgs.")
    
    # In subthreshold, Ids is exponential with Vgs. A +100mV offset on Vs means Vgs changes by -100mV.
    # The ratio of new current to old current is exp(delta_Vgs / (n*Vt)).
    delta_vgs_D = -sensor_offset
    exponent_D = delta_vgs_D / n_V_T
    current_ratio_D = math.exp(exponent_D)
    # Resistance is inversely proportional to current.
    resistance_change_factor_D = 1 / current_ratio_D

    print(f"For a {sensor_offset*1000:.0f}mV offset, the change in Vgs is {delta_vgs_D*1000:.0f}mV.")
    print(f"The current ratio is calculated as: exp({delta_vgs_D:.3f}V / {n_V_T:.3f}V) = exp({exponent_D:.3f})")
    print(f"Resulting Current Ratio: {current_ratio_D:.4f}")
    print(f"This means the resistance changes by a factor of {resistance_change_factor_D:.2f}, which is unacceptable.")
    print("Conclusion: Any design without bootstrapping (like D) fails the offset tolerance requirement.\n")

    # --- Step 2: Analyze Stability (Comparing A and E) ---
    print("--- Step 2: Impact of Switch Charge Injection on Stability ---")
    print("When switching from 'reset' to 'operate', charge (delta_Q) is injected onto the floating gate capacitor (C).")
    print("This creates a voltage error: delta_V = delta_Q / C, which destabilizes the resistance.\n")
    
    # For strategy A
    voltage_error_A = charge_injection_q / cap_A
    exponent_A = voltage_error_A / n_V_T
    # The resistance change factor is 1/current_change_factor
    resistance_change_factor_A = 1 / math.exp(exponent_A) 

    print("Strategy A uses a small 1pF capacitor for a quick reset.")
    print(f"With an assumed {charge_injection_q*1e15:.0f}fC charge injection:")
    print(f"Voltage Error = {charge_injection_q:.1e}C / {cap_A:.1e}F = {voltage_error_A*1000:.1f}mV.")
    print(f"This error voltage causes a resistance change calculated by: exp(-({voltage_error_A:.4f}V / {n_V_T:.3f}V)) = exp({-exponent_A:.3f})")
    print(f"Resulting Resistance Change Factor: {resistance_change_factor_A:.3f} (a {100*(1-resistance_change_factor_A):.1f}% error from just one event).")
    print("This high sensitivity to non-idealities makes Strategy A unstable.\n")

    # For strategy E
    voltage_error_E = charge_injection_q / cap_E
    exponent_E = voltage_error_E / n_V_T
    resistance_change_factor_E = 1 / math.exp(exponent_E)

    print("Strategy E's split-gate allows a fast reset even with a larger 10pF capacitor.")
    print(f"With the same {charge_injection_q*1e15:.0f}fC charge injection:")
    print(f"Voltage Error = {charge_injection_q:.1e}C / {cap_E:.1e}F = {voltage_error_E*1000:.1f}mV.")
    print(f"This error voltage causes a resistance change calculated by: exp(-({voltage_error_E:.4f}V / {n_V_T:.3f}V)) = exp({-exponent_E:.3f})")
    print(f"Resulting Resistance Change Factor: {resistance_change_factor_E:.3f} (a {100*(1-resistance_change_factor_E):.1f}% error).")
    print("This is far more robust and stable than Strategy A.\n")

    # --- Final Conclusion ---
    print("--- Overall Conclusion ---")
    print("Strategy E provides the most effective balance:")
    print("1. Rapid Reset: Achieved by using both gate halves to create a low-resistance path.")
    print("2. Offset Tolerance: Retains the essential bootstrapping mechanism, which is robust to DC offsets.")
    print("3. Stability: Can use a larger capacitor, making it 10x less sensitive to charge injection errors than Strategy A.")
    print("4. Mismatch Compensation: The static half of the gate provides a unique, built-in mechanism to trim for process variations and precisely set the bias point.")

analyze_design_strategies()

print("<<<E>>>")