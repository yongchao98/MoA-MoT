import sys

def calculate_sum_approximation():
    """
    Calculates and prints the approximation for the sum for a given n.
    The value of n is taken from the command line, with a default of 100.
    """
    # The user can specify the value of n as a command-line argument.
    if len(sys.argv) > 1:
        try:
            n = int(sys.argv[1])
            if n <= 0:
                print("Error: n must be a positive integer.", file=sys.stderr)
                sys.exit(1)
        except ValueError:
            print("Error: Invalid integer provided for n.", file=sys.stderr)
            sys.exit(1)
    else:
        # Default value for n if not provided
        n = 100

    # The derived formula for the sum with an error of O(n^{-2}) is:
    # S_n ≈ n^2/2 + 1/120 + 1/(252*n)

    # Calculate the value of each term in the formula
    term1 = n**2 / 2.0
    term2 = 1.0 / 120.0
    term3 = 1.0 / (252.0 * n)

    # Calculate the final approximation
    approximation = term1 + term2 + term3

    # As requested, we output each number in the final equation.
    # This is interpreted as showing the formula with n substituted,
    # and the value of each term that is summed.
    print(f"For n = {n}, the approximation formula is:")
    print(f"S_n ≈ ({n}**2 / 2) + (1 / 120) + (1 / (252 * {n}))")
    print("\nCalculating each term:")
    print(f"Term 1 (n^2/2): {term1}")
    print(f"Term 2 (1/120): {term2}")
    print(f"Term 3 (1/(252*n)): {term3}")
    print("\nFinal equation with calculated values:")
    print(f"{term1} + {term2} + {term3} = {approximation}")
    print(f"\nFinal approximated sum: {approximation}")

if __name__ == "__main__":
    calculate_sum_approximation()