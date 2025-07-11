import sympy

def calculate_charges_formulas():
    """
    This function defines and prints the symbolic formulas for the charges as given in Option B.
    """
    # Define the symbols
    V, epsilon, pi, L, a, b = sympy.symbols('V epsilon pi L a b', real=True, positive=True)

    # Formulas from Option B
    # Total volume charge
    q_v_expr = -4 * V * epsilon * pi * L / (1 - a**2 / b**2)

    # Total surface charge on the inner electrode
    q_s_a_expr = 2 * pi * L * V * epsilon / (1 - a**2 / b**2)

    # Total surface charge on the outer electrode
    # The original formula can be simplified for clarity in printing
    q_s_b_expr = -4 * pi * L * V * epsilon * a**2 / (b**2 * (1 - a**2 / b**2))
    # Simplified form: -4*pi*L*V*epsilon*a**2 / (b**2 - a**2)

    # Print the formulas
    print("Based on the selected answer choice B, the formulas are:")
    print("-" * 50)
    print("Total volume charge (q_v):")
    sympy.pprint(q_v_expr)
    print("\n" + "-" * 50)
    print("Total surface charge on inner electrode q_s(r=a):")
    sympy.pprint(q_s_a_expr)
    print("\n" + "-" * 50)
    print("Total surface charge on outer electrode q_s(r=b):")
    sympy.pprint(q_s_b_expr)
    print("\n" + "-" * 50)


# Execute the function
calculate_charges_formulas()