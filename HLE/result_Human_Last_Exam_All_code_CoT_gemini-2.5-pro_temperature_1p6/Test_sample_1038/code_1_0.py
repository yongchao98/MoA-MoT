import math

def analyze_body_biasing():
    """
    Analyzes the effect of increasing threshold voltage (Vt) on subthreshold current.
    This demonstrates the core principle of option C.
    """
    
    # --- Assumptions ---
    # n: subthreshold slope factor
    # U_T: thermal voltage at room temperature (approx. 26mV)
    # Vt_nominal: The transistor's standard threshold voltage
    # Vgs_bias: The gate-source voltage chosen to bias the transistor in subthreshold
    # Vt_increase: The increase in Vt from the on-chip body-bias generator
    n = 1.5
    U_T = 0.026  # Volts
    Vt_nominal = 0.45  # Volts
    Vgs_bias = 0.40  # Volts (a bias point just below the nominal Vt)
    Vt_increase = 0.20 # Volts (a plausible increase from body biasing)

    print("--- Analysis of Strategy C: Using Body Bias to Increase Vt ---\n")
    print(f"Goal: Reduce subthreshold current to decrease leakage and increase resistance.")
    print(f"The subthreshold current is proportional to: exp((Vgs - Vt) / (n * U_T))\n")

    print("--- Case 1: Nominal Operation (No Body Bias) ---")
    
    # Calculate Vt during nominal operation (no change)
    Vt_case1 = Vt_nominal
    
    # Calculate the denominator of the exponent term
    n_Ut_product = n * U_T

    # Calculate the exponent for the nominal case
    exponent_case1 = (Vgs_bias - Vt_case1) / n_Ut_product
    
    # Calculate a relative measure of the subthreshold current
    relative_current_case1 = math.exp(exponent_case1)

    print(f"Nominal Vt = {Vt_case1:.2f} V")
    print(f"Chosen Vgs Bias = {Vgs_bias:.2f} V")
    print(f"Subthreshold slope factor 'n' = {n}")
    print(f"Thermal Voltage 'U_T' = {U_T} V\n")
    
    print("Calculation of the relative subthreshold current:")
    print(f"  Proportional to: exp( (Vgs - Vt) / (n * U_T) )")
    print(f"  Proportional to: exp( ({Vgs_bias:.2f} - {Vt_case1:.2f}) / ({n} * {U_T}) )")
    print(f"  Proportional to: exp( {Vgs_bias - Vt_case1:.2f} / {n_Ut_product:.3f} )")
    print(f"  Proportional to: exp( {exponent_case1:.3f} )")
    print(f"  Relative Current = {relative_current_case1:.4e}\n")
    
    print("--- Case 2: Operate Phase (Body Bias ON) ---")
    
    # Calculate Vt with body biasing
    Vt_case2 = Vt_nominal + Vt_increase

    # Calculate the exponent for the biased case
    exponent_case2 = (Vgs_bias - Vt_case2) / n_Ut_product

    # Calculate the relative subthreshold current
    relative_current_case2 = math.exp(exponent_case2)
    
    print(f"Body-biased Vt = Nominal Vt + Increase = {Vt_nominal:.2f} V + {Vt_increase:.2f} V = {Vt_case2:.2f} V")
    print(f"Chosen Vgs Bias = {Vgs_bias:.2f} V (remains the same)\n")
    
    print("Calculation of the new relative subthreshold current:")
    print(f"  Proportional to: exp( (Vgs - Vt) / (n * U_T) )")
    print(f"  Proportional to: exp( ({Vgs_bias:.2f} - {Vt_case2:.2f}) / ({n} * {U_T}) )")
    print(f"  Proportional to: exp( {Vgs_bias - Vt_case2:.2f} / {n_Ut_product:.3f} )")
    print(f"  Proportional to: exp( {exponent_case2:.3f} )")
    print(f"  Relative Current = {relative_current_case2:.4e}\n")

    # --- Conclusion ---
    reduction_factor = relative_current_case1 / relative_current_case2
    print("--- Conclusion ---")
    print(f"By increasing Vt from {Vt_case1:.2f}V to {Vt_case2:.2f}V, the subthreshold current is reduced by a factor of {reduction_factor:.1f}.")
    print("This directly translates to lower gate-cap leakage and a higher, more stable resistance during operation, effectively solving the design challenge.")

if __name__ == '__main__':
    analyze_body_biasing()