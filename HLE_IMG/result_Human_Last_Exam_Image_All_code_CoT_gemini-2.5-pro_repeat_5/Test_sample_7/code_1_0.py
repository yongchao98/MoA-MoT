import sympy

def solve_electrochemistry_formula():
    """
    This function derives and prints the formula for the second voltage plateau
    of a graphite anode based on the provided physical description.
    """

    # Given values for the two voltage plateaus
    V1_val = 0.09  # Voltage for the first plateau (Stage 2 -> 1)
    V2_val = 0.13  # Voltage for the second plateau (Stage 3 -> 2)

    # Define symbols for the formula
    V1, V2, mu1, mu2, mu_Li, e = sympy.symbols('V_1 V_2 μ_1 μ_2 μ_Li e')

    # The fundamental equations relating voltage (V) to chemical potential (μ)
    # V = -(μ - μ_Li) / e  which can be written as V*e = μ_Li - μ
    eq1 = sympy.Eq(V1 * e, mu_Li - mu1)
    eq2 = sympy.Eq(V2 * e, mu_Li - mu2)

    # We want to find a formula for V2 in terms of mu1, mu2, and e.
    # We can express V2 in terms of V1 and the chemical potentials by eliminating mu_Li.
    # From eq1, we can solve for mu_Li: mu_Li = V1*e + mu1
    # Substitute this into eq2: V2*e = (V1*e + mu1) - mu2
    # V2 = V1 + (mu1 - mu2)/e
    # This formula expresses the second plateau voltage (V2) in terms of the first plateau voltage (V1)
    # and the two chemical potentials.

    # The problem asks for a formula that best approximates the second plateau.
    # We use the derived relationship and substitute the given value for V1.
    final_formula_rhs = f"{V1_val} + (μ_1 - μ_2)/e"

    print("The simple formula that best approximates the second plateau (V_2) is:")
    print(f"V_2 = V_1 + (μ_1 - μ_2)/e")
    print("\nSubstituting the given numerical value for the first plateau (V_1):")
    # The final equation includes all the numbers as requested.
    # The value of the second plateau V_2 is on the left.
    # The formula on the right uses the value of the first plateau V_1.
    print(f"{V2_val} = {final_formula_rhs}")
    
    # We can verify that this is consistent:
    # V2 - V1 = (μ_1 - μ_2)/e
    # 0.13 - 0.09 = 0.04, so (μ_1 - μ_2)/e = 0.04 V.
    # The formula holds true: 0.13 = 0.09 + 0.04.

solve_electrochemistry_formula()