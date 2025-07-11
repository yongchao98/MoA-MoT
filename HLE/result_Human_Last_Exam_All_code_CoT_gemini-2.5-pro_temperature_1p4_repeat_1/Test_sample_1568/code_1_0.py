def solve_and_print_formula():
    """
    This function constructs and prints the closed-form expression
    for the infinite product P = Product_{n=3 to infinity} (1 - z^3/n^3).
    """

    # Define the numbers that appear in the product definition and final formula.
    start_index = 3
    power = 3
    n1 = 1
    n2 = 2
    root_of_unity_numerator = 2
    root_of_unity_denominator = 3

    # Construct the left-hand side (LHS) of the equation
    lhs = f"Product_{{n={start_index} to infinity}} (1 - z^{power}/n^{power})"

    # Construct the right-hand side (RHS) of the equation
    term1_simple = f"(1 - z^{power})"
    term2 = f"(1 - z^{power}/{n2**power})"
    gamma1 = f"Gamma(1 - z)"
    gamma2 = f"Gamma(1 - z*exp(i*{root_of_unity_numerator}*pi/{root_of_unity_denominator}))"
    gamma3 = f"Gamma(1 - z*exp(-i*{root_of_unity_numerator}*pi/{root_of_unity_denominator}))"

    rhs_numerator = "1"
    rhs_denominator = f"{term1_simple} * {term2} * {gamma1} * {gamma2} * {gamma3}"

    # Print the final equation
    print(f"The value of the infinite product is given by the equation:")
    print()
    print(f"{lhs} = \n")
    print(f"    {rhs_numerator}")
    print(f"    ----------------------------------------------------------------------------------------------------")
    print(f"    {rhs_denominator}")
    print()
    print("Where 'Gamma' is the Gamma function, 'exp' is the exponential function,")
    print("'i' is the imaginary unit, and 'pi' is the mathematical constant pi.")
    print()

    # Per instructions, outputting each number in the final equation.
    print("The numbers that appear in this final equation are:")
    # Numbers are from the LHS
    print(f"Start index of product: {start_index}")
    print(f"Power of z and n: {power}")
    # Numbers from the RHS
    print(f"Numerator: {n1}")
    print(f"Exponent in the (1 - z^3) term: {power}")
    print(f"Base in the (1 - z^3/8) term: {n2}")
    print(f"Resulting denominator ({n2}^{power}): {n2**power}")
    print(f"Constant term in Gamma argument '1 - z': 1")
    print(f"Numerator for pi in exp argument: {root_of_unity_numerator}")
    print(f"Denominator for pi in exp argument: {root_of_unity_denominator}")

# Execute the function to print the solution.
solve_and_print_formula()