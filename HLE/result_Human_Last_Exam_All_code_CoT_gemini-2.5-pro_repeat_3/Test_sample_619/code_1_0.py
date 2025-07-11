import math

def find_minimum_sum(n):
    """
    Calculates the minimum value of the sum of cardinalities of n sets S_i
    satisfying |S_i triangle S_j| = |i-j|.

    Args:
        n: The number of sets.

    Returns:
        The minimum possible value of sum(|S_i|).
    """
    if not isinstance(n, int) or n < 1:
        raise ValueError("n must be a positive integer.")
    
    # The derived formula for the minimum sum is floor(n^2 / 4) + n.
    # For n=2k, this is k^2 + 2k.
    # For n=2k+1, this is k^2 + 3k + 1.
    
    term1_val = math.floor(n**2 / 4)
    term2_val = n
    min_sum = term1_val + term2_val
    
    return min_sum

# You can replace this value with any integer n >= 1 to test.
n = 7

# Calculate the result
min_value = find_minimum_sum(n)

# Print the final equation and result
print(f"For n = {n}:")
print(f"The minimum value of the sum is floor({n}^2 / 4) + {n} = floor({n**2}/4) + {n} = {math.floor(n**2/4)} + {n} = {min_value}")

# Comparing this result with the answer choices for n=7:
# A. 12
# B. 14
# C. 50
# D. 49
# E. 13
# F. None of the above
# The calculated value 19 is not among choices A-E.