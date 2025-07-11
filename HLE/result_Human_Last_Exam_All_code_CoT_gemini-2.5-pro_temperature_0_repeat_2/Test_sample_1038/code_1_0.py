import math

def demonstrate_body_effect():
    """
    This script demonstrates the principle of Option C by calculating how
    reverse body-biasing increases a transistor's threshold voltage (Vt)
    and reduces its subthreshold current, thus increasing its effective resistance.
    """

    # --- Assumed Transistor Parameters ---
    Vt0 = 0.45  # Intrinsic threshold voltage in Volts (from prompt)
    gamma = 0.4  # Body effect coefficient in V^0.5 (typical value)
    phi_f = 0.35 # Fermi potential in Volts (typical value, so 2*phi_f is ~0.7V)
    n = 1.5      # Subthreshold slope factor (from prompt)
    U_T = 0.026  # Thermal voltage in Volts (~room temperature)
    I0 = 1e-18   # Process-dependent constant for subthreshold current in Amps (example value)
    Vgs = 0.5    # An example gate-source voltage in Volts

    print("--- Transistor Parameters ---")
    print(f"Intrinsic Threshold Voltage (Vt0): {Vt0} V")
    print(f"Body Effect Coefficient (gamma): {gamma} V^0.5")
    print(f"Surface Potential (2*phi_f): {2*phi_f} V")
    print(f"Subthreshold Slope Factor (n): {n}")
    print(f"Assumed Gate-Source Voltage (Vgs): {Vgs} V\n")

    # --- Case 1: No Body Bias (Standard Operation) ---
    Vbs1 = 0.0  # Body-source voltage in Volts
    
    # Calculate Vt using the body effect equation: Vt = Vt0 + gamma * (sqrt(2*phi_f - Vbs) - sqrt(2*phi_f))
    sqrt_term1 = math.sqrt(2 * phi_f - Vbs1)
    Vt1 = Vt0 + gamma * (sqrt_term1 - math.sqrt(2 * phi_f))
    
    # Calculate subthreshold current: Id = I0 * exp((Vgs - Vt) / (n * U_T))
    exponent1 = (Vgs - Vt1) / (n * U_T)
    Id1 = I0 * math.exp(exponent1)
    
    print("--- Case 1: No Body Bias (Vbs = 0V) ---")
    print(f"Equation for Vt: {Vt0} + {gamma} * (sqrt({2*phi_f} - {Vbs1}) - sqrt({2*phi_f}))")
    print(f"Calculated Threshold Voltage (Vt1): {Vt1:.4f} V")
    print(f"Equation for Id exponent: ({Vgs} - {Vt1:.4f}) / ({n} * {U_T}) = {exponent1:.4f}")
    print(f"Resulting Subthreshold Current (Id1): {Id1:.4e} A\n")

    # --- Case 2: With Reverse Body Bias (Option C) ---
    Vbs2 = -0.3 # Reverse body bias in Volts (as suggested in Option C for NMOS)

    # Recalculate Vt with the new Vbs
    sqrt_term2 = math.sqrt(2 * phi_f - Vbs2)
    Vt2 = Vt0 + gamma * (sqrt_term2 - math.sqrt(2 * phi_f))

    # Recalculate subthreshold current
    exponent2 = (Vgs - Vt2) / (n * U_T)
    Id2 = I0 * math.exp(exponent2)

    print("--- Case 2: With Reverse Body Bias (Vbs = -0.3V) ---")
    print(f"Equation for Vt: {Vt0} + {gamma} * (sqrt({2*phi_f} - ({Vbs2})) - sqrt({2*phi_f}))")
    print(f"Calculated Threshold Voltage (Vt2): {Vt2:.4f} V")
    print(f"Equation for Id exponent: ({Vgs} - {Vt2:.4f}) / ({n} * {U_T}) = {exponent2:.4f}")
    print(f"Resulting Subthreshold Current (Id2): {Id2:.4e} A\n")

    # --- Comparison and Conclusion ---
    resistance_increase_factor = Id1 / Id2
    print("--- Conclusion ---")
    print(f"Applying a {Vbs2}V body bias increased Vt by {Vt2 - Vt1:.4f} V.")
    print(f"This reduced the subthreshold current by a factor of {resistance_increase_factor:.2f}.")
    print("A lower current for the same Vgs means a higher effective resistance, demonstrating the effectiveness of Option C.")

if __name__ == '__main__':
    demonstrate_body_effect()