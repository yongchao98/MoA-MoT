import sympy as sp

def display_force_formula():
    """
    This function displays the components of the derived force per unit area formula.
    """
    # Define symbolic variables
    mu_0, mu, K_0, a, y, d = sp.symbols('mu_0 mu K_0 a y d')
    
    # Numerator of the main fraction
    numerator_main = K_0**2 * sp.sin(a*y)**2
    
    # Denominator of the main fraction
    # First term inside the bracket
    denom_term1 = sp.cosh(a*d)
    # Second term inside the bracket (fraction)
    denom_term2_num = mu_0
    denom_term2_den = mu
    denom_term2 = (denom_term2_num / denom_term2_den) * sp.sinh(a*d)
    # Full denominator squared
    denominator_main = (denom_term1 + denom_term2)**2

    # Prefactor fraction
    prefactor_num = mu_0
    prefactor_den = 2

    # Construct the final formula string for printing
    print("The derived force per unit area is:")
    print("f/area = (prefactor) * (numerator) / (denominator) * i_x")
    print("\nWhere the components are:")
    
    print(f"\nPrefactor Numerator: {prefactor_num}")
    # Printing each "number" or symbol as requested
    print(f"The number in the prefactor denominator is: 2")

    print(f"\nMain Fraction Numerator: {numerator_main}")

    print(f"\nMain Fraction Denominator: ({sp.pretty(denom_term1)} + ({sp.pretty(denom_term2_num)}/{sp.pretty(denom_term2_den)})*{sp.pretty(sp.sinh(a*d))})**2")

    print("\nSo the full expression is:")
    full_expression = (prefactor_num / prefactor_den) * (numerator_main / denominator_main)
    # Using sp.pretty for a more readable output
    print(sp.pretty(full_expression) + " * i_x")


display_force_formula()