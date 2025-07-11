import math

def combinations(n, k):
    """Calculates the number of combinations of k items from a set of n."""
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

# The number of special points in the metric continuum.
num_points = 5

# Each subcontinuum in the decomposition can contain at most 2 of the special points.
# To maximize n, we should use subcontinua that correspond to pairs of points.
# The maximum number n is the number of ways to choose 2 points from the 5 available points.
max_size_of_subset = 2

# Calculate the number of combinations.
n = combinations(num_points, max_size_of_subset)

print(f"The problem is to find the largest number of subcontinua, n, that can form an irreducible cover of the space X.")
print(f"The conditions on X imply that each subcontinuum can contain at most 2 of the 5 special points.")
print(f"The maximum number n is therefore the number of ways to choose pairs of points from the 5 points.")
print(f"This is a combinatorial calculation: C(5, 2).")
print(f"C(5, 2) = 5! / (2! * (5-2)!)")
print(f"5! = {math.factorial(5)}")
print(f"2! = {math.factorial(2)}")
print(f"3! = {math.factorial(3)}")
print(f"n = {math.factorial(5)} / ({math.factorial(2)} * {math.factorial(3)}) = {n}")
print(f"\nThe largest possible number n is {n}.")
