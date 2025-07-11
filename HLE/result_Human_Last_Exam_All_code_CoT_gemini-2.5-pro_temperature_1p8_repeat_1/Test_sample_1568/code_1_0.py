def solve_infinite_product():
    """
    This function constructs and prints the symbolic expression for the infinite product
    prod_{n=3 to inf} (1 - z^3 / n^3).
    """
    # The numbers that appear in the simplified final equation are 8, 3, 2, and 3.
    # We will construct the string representation of the formula.
    numerator = 8
    gamma_arg_const = 3
    omega_num = 2
    omega_den = 3

    # Symbolic variables
    z = "z"
    Gamma = "Gamma"
    pi = "pi"
    i = "i"

    # String for omega and omega squared
    w = f"exp({omega_num}*{pi}*{i}/{omega_den})"
    w_squared = f"exp({2*omega_num}*{pi}*{i}/{omega_den})"

    # Building the denominator string part by part
    term1 = f"{Gamma}({gamma_arg_const} - {z})"
    term2 = f"{Gamma}({gamma_arg_const} - {z}*{w})"
    term3 = f"{Gamma}({gamma_arg_const} - {z}*{w_squared})"

    denominator = f"{term1} * {term2} * {term3}"

    # Final equation as a string
    final_equation = f"{numerator} / ({denominator})"

    # Print the result.
    # The output shows the final equation with all its numbers.
    print(final_equation)

solve_infinite_product()