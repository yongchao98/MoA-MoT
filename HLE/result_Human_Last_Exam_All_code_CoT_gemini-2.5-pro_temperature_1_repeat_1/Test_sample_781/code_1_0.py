import math

# The number of special points in the continuum X.
num_points = 5

# As explained in the reasoning, any proper subcontinuum in the decomposition
# can contain at most 2 of the special points.
points_per_subcontinuum = 2

# The largest number n is the number of ways to choose pairs of points from the
# set of special points. This is calculated using the combination formula C(n, k).
result = math.comb(num_points, points_per_subcontinuum)

# For displaying the full equation, we calculate the factorials.
fact_total = math.factorial(num_points)
fact_chosen = math.factorial(points_per_subcontinuum)
fact_remainder = math.factorial(num_points - points_per_subcontinuum)

print("The largest number n is the number of pairs that can be chosen from the 5 special points.")
print("This is calculated by the combination formula: C(m, k) = m! / (k! * (m-k)!).")
print(f"In this problem, m = {num_points} and k = {points_per_subcontinuum}.")
print("\nThe final equation with the numbers is:")
print(f"{result} = {fact_total} / ({fact_chosen} * {fact_remainder})")