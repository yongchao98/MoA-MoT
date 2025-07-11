import math

def calculate_combinations(n, k):
    """Calculates the number of combinations of choosing k items from a set of n."""
    return math.comb(n, k)

# The problem states there are 5 distinct points.
num_points = 5

# Each subcontinuum in the decomposition can contain at most 2 of these points.
# The largest number of such continua corresponds to the number of distinct pairs
# of points that can be formed from the set of 5 points.
num_chosen = 2

# Calculate the number of combinations.
result = calculate_combinations(num_points, num_chosen)

# Output the explanation and the result.
print(f"The problem asks for the largest number n such that X can be decomposed into n subcontinua A_i,")
print(f"where each A_i is a proper subcontinuum and thus contains at most 2 of the 5 special points.")
print(f"This number is maximized when each subcontinuum corresponds to a unique pair of the special points.")
print(f"The number of unique pairs from 5 points is the combination C(5, 2).")
print(f"C({num_points}, {num_chosen}) = {result}")
