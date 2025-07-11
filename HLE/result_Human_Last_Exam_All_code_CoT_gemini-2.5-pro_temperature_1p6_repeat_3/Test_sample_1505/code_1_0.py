import sys

def solve_sum_approximation():
    """
    Calculates the approximation for the sum for a given integer n.
    The user should set the value of n below.
    """
    # Set the value for the positive integer n.
    # The user can change this value. For example, n = 100.
    try:
        if len(sys.argv) > 1:
            n = int(sys.argv[1])
        else:
            # Fallback to a default value if no command-line argument is given.
            n = 10
            print(f"No value for n provided. Using default value n = {n}.")
            print("You can also run this script from the command line, e.g., 'python your_script.py 100'")

    except ValueError:
        print("Invalid input. Please provide a positive integer for n.", file=sys.stderr)
        sys.exit(1)

    # The approximation formula for the sum is n^2/2 + 1/120 + 1/(252*n).
    # We calculate each part of this formula, which are the numbers in the equation.
    term1 = float(n * n) / 2.0
    term2 = 1.0 / 120.0
    term3 = 1.0 / (252.0 * float(n))

    # The result of the equation is the sum of the terms.
    approximation = term1 + term2 + term3

    # As requested, print each number in the final equation.
    # The "equation" is: term1 + term2 + term3 = approximation.
    # The numbers are the values of term1, term2, term3, and the final result.
    print(f"For n = {n}, the terms of the approximation formula are:")
    print(f"Term 1 (n^2 / 2): {term1}")
    print(f"Term 2 (1 / 120): {term2}")
    print(f"Term 3 (1 / (252 * n)): {term3}")
    print(f"The final approximation is: {approximation}")

if __name__ == "__main__":
    solve_sum_approximation()