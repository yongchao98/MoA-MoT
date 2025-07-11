import math

def analyze_pseudo_resistor_strategy():
    """
    Analyzes the design constraints for a bootstrapped pseudo-resistor
    and determines the most effective strategy among the given options.
    """

    # --- Problem Parameters ---
    VDD = 1.2  # Volts, supply voltage
    Vt0 = 0.45  # Volts, baseline threshold voltage
    V_offset = 0.1  # Volts, sensor offset range
    Vcm = VDD / 2  # Volts, common-mode voltage
    
    # Body effect parameters (typical for CMOS)
    gamma = 0.5  # V^0.5, body effect coefficient
    phi_f2 = 0.7  # Volts, 2 * Fermi potential

    print("--- Analysis of Pseudo-Resistor Design Strategies ---\n")
    print(f"Circuit Constraints: VDD={VDD}V, Vt={Vt0}V, Sensor Offset=+/={V_offset}V\n")

    # --- 1. Baseline Analysis (Without special strategy) ---
    print("1. Baseline Scenario Analysis (No Body Binsing)")
    Vs_max = Vcm + V_offset
    # To maintain subthreshold operation, Vgs must be less than Vt. Vg - Vs < Vt
    # This means the gate voltage Vg must be less than Vs + Vt.
    Vg_limit_baseline = Vs_max + Vt0
    
    print(f"Worst-case source voltage (Vs_max) due to offset: {Vcm:.1f}V + {V_offset:.1f}V = {Vs_max:.1f}V")
    print("To keep the NMOS in subthreshold, Vg < Vs_max + Vt0")
    print(f"Required Vg < {Vs_max:.1f}V + {Vt0:.2f}V = {Vg_limit_baseline:.2f}V")
    print(f"This leaves a very small margin ({VDD - Vg_limit_baseline:.2f}V) below the {VDD}V supply rail, making the bias point sensitive and difficult to maintain.\n")

    # --- 2. Strategy C Analysis (With On-Chip Body Bias) ---
    print("2. Strategy C Analysis (Using Switched Reverse Body Bias)")
    # During the 'operate' phase, apply reverse body bias to the NMOS.
    # We connect the NMOS body to ground (0V). The source is at Vs_max.
    # The body-to-source voltage (Vbs) is negative.
    Vbs = 0 - Vs_max
    
    # Calculate the new, increased threshold voltage (Vt_eff) using the body effect formula.
    # Vt_eff = Vt0 + gamma * (sqrt(2*phi_f - Vbs) - sqrt(2*phi_f))
    sqrt_term_no_bias = math.sqrt(phi_f2)
    sqrt_term_with_bias = math.sqrt(phi_f2 - Vbs)
    Vt_increase = gamma * (sqrt_term_with_bias - sqrt_term_no_bias)
    Vt_eff = Vt0 + Vt_increase
    
    print(f"At Vs_max = {Vs_max:.1f}V, the reverse body-bias Vbs = 0V - {Vs_max:.1f}V = {Vbs:.1f}V")
    print(f"This increases the threshold voltage. The new effective Vt (Vt_eff) is calculated:")
    print(f"Vt_eff = Vt0 + γ * (sqrt(2φf - Vbs) - sqrt(2φf))")
    print(f"Vt_eff = {Vt0:.2f} + {gamma:.1f} * (sqrt({phi_f2:.1f} - ({Vbs:.1f})) - sqrt({phi_f2:.1f}))")
    print(f"Vt_eff = {Vt0:.2f} + {gamma:.1f} * ({sqrt_term_with_bias:.2f} - {sqrt_term_no_bias:.2f})")
    print(f"Vt_eff = {Vt0:.2f} + {Vt_increase:.2f} = {Vt_eff:.2f}V\n")
    
    # Recalculate the required gate voltage limit with the new Vt_eff.
    Vg_limit_strategy_C = Vs_max + Vt_eff
    
    print("With the higher Vt, the condition for subthreshold operation (Vg < Vs_max + Vt_eff) becomes:")
    
    vs_max_val = Vs_max
    vt_eff_val = Vt_eff
    result_val = Vg_limit_strategy_C
    print(f"Required Vg < {vs_max_val:.1f}V + {vt_eff_val:.2f}V = {result_val:.2f}V")
    
    print(f"Since {Vg_limit_strategy_C:.2f}V is well above the {VDD}V supply, any gate voltage is guaranteed to keep the device in subthreshold.")
    print("This dramatically increases the operating margin and reduces subthreshold leakage current.\n")

    # --- 3. Checking the Reset Phase for Strategy C ---
    print("3. Reset Phase Check for Strategy C")
    print("During reset, body bias is removed (Vbs=0), so Vt reverts to Vt0 = 0.45V.")
    # For a fast reset, the transistor must be in strong inversion (Vgs > Vt0).
    # Assume the reset switch pulls the gate to VDD.
    Vg_reset = VDD
    Vgs_reset = Vg_reset - Vs_max
    print(f"With Vg pulled to {VDD}V and worst-case Vs at {Vs_max:.1f}V, Vgs = {Vg_reset:.1f}V - {Vs_max:.1f}V = {Vgs_reset:.1f}V.")
    
    is_strong_inversion = Vgs_reset > Vt0
    print(f"Is Vgs > Vt0?  {Vgs_reset:.1f}V > {Vt0:.2f}V? -> {is_strong_inversion}")
    print("Strong inversion is achieved, ensuring a fast reset is still possible.\n")

    # --- 4. Conclusion ---
    print("Conclusion: Strategy C provides the best balance. It uses switched body-biasing to increase Vt during the high-impedance 'operate' phase, which solves the headroom and leakage problems. It then reverts to the lower, nominal Vt during the 'reset' phase, allowing for fast settling. The other options introduce significant new problems like leakage, noise, or functional limitations.")

analyze_pseudo_resistor_strategy()