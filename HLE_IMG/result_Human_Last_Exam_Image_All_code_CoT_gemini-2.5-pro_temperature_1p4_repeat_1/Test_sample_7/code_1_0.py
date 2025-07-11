import sympy

def solve_and_print_equation():
    """
    This function derives and prints the formula for the second voltage plateau
    based on the provided physical parameters.
    """
    # Given empirical values for the two voltage plateaus
    V1 = 0.09  # Voltage of the first plateau (Stage 2 -> 1) in Volts
    V2 = 0.13  # Voltage of the second plateau (Stage 3 -> 2) in Volts

    # Define symbolic variables for the chemical potentials and charge
    mu1, mu2, e = sympy.symbols("μ_1 μ_2 e")
    
    # The voltage difference between the two plateaus is related to the difference
    # in chemical potentials.
    # V2 - V1 = ((μ_Li - μ_2)/e) - ((μ_Li - μ_1)/e) = (μ_1 - μ_2)/e
    voltage_difference_term = V2 - V1

    # The formula for the second plateau V2 can be expressed as:
    # V2 = V1 + (V2 - V1)
    # Substituting the physical meaning of the difference:
    # V2 = V1 + (μ_1 - μ_2)/e
    
    # Print the explanation and the final formula with numerical values.
    print("The relationship between the two plateau voltages V1 and V2 and the chemical potentials μ1 and μ2 is:")
    print("V2 = V1 + (μ_1 - μ_2) / e")
    print("\nPlugging in the numerical values given in the problem:")
    print(f"{V2:.2f} V = {V1:.2f} V + ({voltage_difference_term:.2f} V)")
    print("\nWhere the term (μ_1 - μ_2) / e, which represents the difference in the plateau voltages, is calculated as:")
    print(f"{V2:.2f} V - {V1:.2f} V = {voltage_difference_term:.2f} V")
    
    # The final symbolic formula is the most direct answer to the question.
    final_formula_str = f"V_1 + (μ_1 - μ_2)/e"
    print(f"\nThe symbolic formula for the second plateau is: {final_formula_str}")


solve_and_print_equation()