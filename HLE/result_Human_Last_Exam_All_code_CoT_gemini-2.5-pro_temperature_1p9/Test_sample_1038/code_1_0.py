import math

def analyze_option_a():
    """
    Analyzes the quantitative trade-offs of the pseudo-resistor design
    described in Option A.
    """
    print("--- Analysis of Design Strategy A ---")
    print("Strategy A: Minimum-length, large-width transistors with a small gate capacitor (1 pF).\n")

    # --- Given Parameters and Assumptions ---
    C_gate = 1e-12  # Gate capacitance = 1 pF
    t_reset_spec = 5e-6  # Reset time specification < 5 µs
    leakage_spec = 0.01  # Gate voltage droop < 1% per second
    Vt = 0.45  # Threshold voltage = 0.45 V
    Vdd = 1.2  # Supply voltage = 1.2 V
    
    # Assume the source is biased mid-rail for the reset calculation
    Vs_reset = 0.6 # V
    Vgs_reset = Vdd - Vs_reset # Vgs during reset is Vdd - Vs
    
    # Transistor parameter assumptions for calculation
    # Let's assume a large W/L ratio for a strong 'on' transistor, e.g., 100
    W_div_L = 100.0
    # Process transconductance parameter for a modern CMOS process
    k_prime = 200e-6  # A/V^2
    
    # --- 1. Reset Time Analysis ---
    print("1. Checking Reset Time (< 5 microseconds)")
    R_on = 1 / (k_prime * W_div_L * (Vgs_reset - Vt))
    print(f"   - On-Resistance (R_on) = 1 / (k' * W/L * (Vgs - Vt))")
    print(f"   - R_on = 1 / ({k_prime:.2e} A/V^2 * {W_div_L} * ({Vgs_reset:.2f} V - {Vt:.2f} V)) = {R_on:.1f} Ohms")

    tau_reset = R_on * C_gate
    print(f"   - RC Time Constant (tau) = R_on * C_gate")
    print(f"   - tau = {R_on:.1f} Ohms * {C_gate:.1e} F = {tau_reset:.3e} seconds")
    
    # Settling time is typically ~5-7 time constants
    t_settle = 7 * tau_reset
    print(f"   - Estimated settling time (7*tau) is {t_settle:.3e} seconds ({t_settle*1e9:.2f} nanoseconds).")

    if t_settle < t_reset_spec:
        print("   - Conclusion: The reset time is extremely fast, well within the < 5 µs specification.\n")
    else:
        print("   - Conclusion: The reset time FAILS to meet the < 5 µs specification.\n")
        
    # --- 2. Leakage Budget Analysis ---
    print("2. Checking Gate Leakage Tolerance (< 1% per second droop)")
    # Assume the gate voltage in operate mode is near Vt to be in subthreshold
    V_gate_op = Vt + 0.1 # Example operating voltage, e.g., 550mV
    
    # The requirement is (dV/V)/dt < 0.01. So I_leak/C_gate < 0.01 * V_gate
    I_leak_max = leakage_spec * C_gate * V_gate_op
    print(f"   - Max leakage current (I_leak) = (Droop Spec) * C_gate * V_gate")
    print(f"   - I_leak_max = {leakage_spec} * {C_gate:.1e} F * {V_gate_op:.2f} V = {I_leak_max:.2e} Amps")
    print(f"   - This means the total leakage current must be less than {I_leak_max*1e15:.2f} femtoamps.")
    print("   - Conclusion: This is an extremely low current, very challenging to achieve with minimum-length transistors, which are prone to higher leakage.\n")

    # --- 3. Charge Injection Analysis ---
    print("3. Checking Bias Error from Charge Injection")
    # Subthreshold slope factor n=1.5, Thermal voltage kT/q = 26mV
    n_factor = 1.5
    V_T = 0.026 # Volts
    
    # Estimate channel charge (Q_ch = C_ox * W * L * (Vgs - Vt))
    # C_ox is ~10fF/um^2. For W/L=100 and L=0.12um, W=12um. Area=1.44um^2.
    C_ox_area = 1.44e-12 # m^2
    C_ox_per_area = 10e-15 / (1e-6*1e-6) # F/m^2 -> 10 fF/um^2
    C_ox_total = C_ox_area * C_ox_per_area

    Q_channel = C_ox_total * (Vgs_reset - Vt)
    # Assume 50% of channel charge is injected onto the gate capacitor
    Q_injected = 0.5 * Q_channel
    
    print(f"   - Channel Charge (Q_ch) = C_ox_total * (Vgs_reset - Vt)")
    print(f"   - Q_ch = {C_ox_total:.2e} F * ({Vgs_reset:.2f} V - {Vt:.2f} V) = {Q_channel:.3e} Coulombs")
    
    delta_V = Q_injected / C_gate
    print(f"   - Injected Charge (Q_inj, ~50% of Q_ch) = {Q_injected:.3e} C")
    print(f"   - Gate Voltage Error (Delta V) = Q_inj / C_gate")
    print(f"   - Delta V = {Q_injected:.3e} C / {C_gate:.1e} F = {delta_V:.4f} V ({delta_V*1000:.2f} mV)")
    
    # Calculate the impact on subthreshold current
    # % Change in current = Delta_V / (n * V_T)
    percent_I_change = (delta_V / (n_factor * V_T)) * 100
    print(f"   - Subthreshold Current I is proportional to exp(Vgs / (n*V_T))")
    print(f"   - The % error in current is approx. Delta_V / (n * V_T)")
    print(f"   - % Error = {delta_V*1000:.2f} mV / ({n_factor} * {V_T*1000:.1f} mV) = {percent_I_change:.1f}%")
    print(f"   - Conclusion: Charge injection creates a significant {percent_I_change:.1f}% error in the pseudo-resistance each time the circuit switches from reset to operate. This makes precise, stable biasing very difficult.\n")

    print("--- Overall Conclusion for Strategy A ---")
    print("Option A provides a very fast reset but suffers from severe drawbacks:")
    print("1. The required gate leakage is extremely low and hard to achieve with min-length devices.")
    print("2. Charge injection causes a significant and unpredictable shift in the operating point.")
    print("This indicates that a more robust strategy is needed to balance all requirements.")

if __name__ == '__main__':
    analyze_option_a()