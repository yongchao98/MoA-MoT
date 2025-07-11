def solve_and_print_product_formula():
    """
    This function constructs and prints the closed-form expression for the
    infinite product prod_{n=3 to inf}(1 - z^3/n^3).
    """

    # Define the components of the final equation to explicitly show all numbers.
    n_start = 3
    power = 3

    # Left-hand side of the equation
    lhs = f"Product_{{n={n_start} to inf}}(1 - z^{power}/n^{power})"

    # Right-hand side components
    term_n1 = f"(1 - z^{power}/{1**power})"
    term_n2 = f"(1 - z^{power}/{2**power})"
    gamma_product = f"Gamma(1 - z) * Gamma(1 - z*w) * Gamma(1 - z*w^2)"

    # Right-hand side of the equation
    rhs = f"1 / ({term_n1} * {term_n2} * {gamma_product})"
    
    # Explanation of symbols used
    explanation = ("where w = e^(2*pi*i/3) is the primitive cube root of unity,"
                   " and Gamma is the Gamma function.")

    # Print the full result
    print(f"{lhs} =")
    print(rhs)
    print()
    print(explanation)

solve_and_print_product_formula()