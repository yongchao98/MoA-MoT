import math

# The number of special points in the continuum X.
num_points = 5

# As deduced from the problem's conditions, each subcontinuum in the decomposition
# can contain at most 2 of these special points. To maximize the number of
# subcontinua, n, we can associate one subcontinuum with each unique pair of points.
# The problem is therefore equivalent to calculating "num_points choose 2".
num_in_pair = 2

# Calculate the numerator for the combination formula: n * (n-1)
numerator = num_points * (num_points - 1)

# Calculate the denominator: k! which is 2! = 2
denominator = num_in_pair

# The largest number n is the result of the combination calculation.
n = numerator / denominator

# Print the final equation with all its components, as requested.
print(f"The largest number n is the number of ways to choose 2 points from 5.")
print(f"n = C({num_points}, {num_in_pair})")
print(f"{int(n)} = ({num_points} * {num_points - 1}) / {denominator}")