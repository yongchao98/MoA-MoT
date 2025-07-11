import math

def calculate_sum_approximation(n):
    """
    Calculates the approximation for the sum S(n) = sum_{k>=0} k^3 * exp(-k^2/n)
    for a given positive integer n.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.")
        return

    # The approximation formula is: n^2/2 + 1/120 + 1/(252*n)
    term1 = n**2 / 2.0
    term2 = 1.0 / 120.0
    term3 = 1.0 / (252.0 * n)

    approximation = term1 + term2 + term3

    # Print the final equation with each number explicitly shown
    print(f"For n = {n}, the sum is approximated by the equation:")
    print(f"S({n}) ≈ {n}^2/2 + 1/120 + 1/(252*{n})")
    print(f"     = {term1} + {term2:.8f} + {term3:.8f}")
    print(f"     = {approximation:.8f}")

if __name__ == '__main__':
    # You can change the value of n here
    n_value = 100
    calculate_sum_approximation(n_value)