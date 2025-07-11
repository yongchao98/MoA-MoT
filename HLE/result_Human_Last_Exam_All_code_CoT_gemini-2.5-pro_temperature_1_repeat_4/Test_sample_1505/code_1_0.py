import math

def solve_sum_approximation():
    """
    Calculates the approximation for the sum S(n) for a given integer n.
    The user should modify the value of n in the line below.
    """
    # Please set the positive integer value for n here.
    n = 10

    if not isinstance(n, int) or n <= 0:
        print("Error: Please provide a positive integer for n.")
        return

    # The formula is S(n) ~ n^2/2 + 1/120 + 1/(252*n)
    term1 = n**2 / 2
    term2 = 1 / 120
    term3 = 1 / (252 * n)
    
    result = term1 + term2 + term3

    print(f"The approximation formula for the sum is: n^2/2 + 1/120 + 1/(252*n)")
    print(f"For the given n = {n}, the terms are:")
    
    # Output each number in the final equation
    print(f"Term 1 (n^2/2): {n}^2 / 2 = {term1}")
    print(f"Term 2 (1/120): 1 / 120 = {term2}")
    print(f"Term 3 (1/(252*n)): 1 / (252 * {n}) = {term3}")
    
    print("\nThe final approximate sum is:")
    print(f"{term1} + {term2} + {term3} = {result}")

solve_sum_approximation()