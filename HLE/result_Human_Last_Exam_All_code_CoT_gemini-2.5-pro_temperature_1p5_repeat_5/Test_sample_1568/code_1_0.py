def solve_infinite_product():
    """
    This function prints the closed-form expression for the given infinite product.
    """

    # The expression for the product from n=1 to infinity involves Gamma functions.
    gamma_product_part = "Gamma(1 - z) * Gamma(1 - z*w) * Gamma(1 - z*w**2)"

    # The original product starts from n=3, so we divide by the terms for n=1 and n=2.
    n1_term_num = 1
    n1_term_power = 3
    n1_term = f"(1 - z**{n1_term_power}/{n1_term_num}**{n1_term_power})"

    n2_term_num = 2
    n2_term_power = 3
    n2_term = f"(1 - z**{n2_term_power}/{n2_term_num}**{n2_term_power})"

    # Construct the final expression string.
    # The numerator is 1.
    numerator = "1"
    
    # The denominator is the product of the terms for n=1, n=2, and the Gamma product.
    denominator = f"{n1_term} * {n2_term} * {gamma_product_part}"

    # The final expression.
    final_expression = f"1 / ({denominator})"

    print(f"The value of the infinite product is:")
    print(final_expression)
    print("\nwhere 'Gamma' is the Gamma function and 'w' is the principal cube root of unity, w = exp(2*pi*i/3).")
    # For clarity, also show the expression with simplified constants.
    simplified_denominator = f"(1 - z**3) * (1 - z**3/8) * {gamma_product_part}"
    simplified_expression = f"1 / ({simplified_denominator})"
    print("\nWith simplified constants, the expression is:")
    print(simplified_expression)

solve_infinite_product()