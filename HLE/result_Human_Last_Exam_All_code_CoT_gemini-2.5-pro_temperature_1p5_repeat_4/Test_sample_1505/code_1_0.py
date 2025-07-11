import math

def calculate_sum_approximation(n):
    """
    Calculates the approximation for the sum S_n for a given positive integer n.
    The approximation formula is S_n ≈ n^2/2 + 1/120 + 1/(252*n).
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.")
        return

    # The coefficients/divisors in the formula
    c1 = 2
    c2 = 120
    c3 = 252

    # Calculate the approximation
    approximation = (n**2) / c1 + 1 / c2 + 1 / (c3 * n)

    # Print the final equation with all numbers, as requested.
    print(f"The approximation formula is: (n^2)/{c1} + 1/{c2} + 1/({c3}*n)")
    print(f"For n = {n}, the result of the equation is:")
    # Using format specifiers to clearly show the equation's structure and result
    print(f"({n}^2)/{c1} + 1/{c2} + 1/({c3}*{n}) = {approximation}")

# You can change the value of n here
n_value = 100
calculate_sum_approximation(n_value)