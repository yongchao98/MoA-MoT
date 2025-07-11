def solve_infinite_product():
    """
    This function constructs and prints the symbolic expression for the infinite product
    from n=3 to infinity of (1 - z^3/n^3).
    """

    # Define the numbers that appear in the equation for clarity and to meet the user's request.
    start_index = 3
    power = 3
    num_1 = 1
    num_2 = 2
    num_2_cubed = num_2**power

    # Build the string representing the final expression.
    # The numerator of the final expression is 1.
    numerator = f"{num_1}"

    # The denominator consists of the terms for n=1 and n=2, and the Gamma function product.
    term_n1 = f"(1 - z^{power})"
    term_n2 = f"(1 - z^{power}/{num_2_cubed})"
    gamma_product_str = f"Gamma({num_1} - z) * Gamma({num_1} - z*w) * Gamma({num_1} - z*w^2)"

    # Combine the parts into the final equation string.
    final_equation = f"{numerator} / ({term_n1} * {term_n2} * {gamma_product_str})"

    # Print the full explanation and the resulting formula.
    print(f"The value of the infinite product Product_{{n={start_index} to inf}}(1 - z^{power}/n^{power}) is given by the following expression:")
    print("\n" + final_equation + "\n")
    print("where:")
    print("  Gamma(x) is the Gamma function.")
    print(f"  w = exp(2*pi*i/{power}), which is the principal cube root of unity.")
    print(f"  w^2 = exp(4*pi*i/{power}), which is the other complex cube root of unity.")

# Execute the function to display the solution.
solve_infinite_product()