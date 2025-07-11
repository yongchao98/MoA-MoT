import math

# Problem parameters
initial_radius = 3  # r_1, in meters

# The visibility constraint, derived from m'(r) > 0, is K < r^2.
# This must hold for the smallest r, which is the initial_radius.
k_upper_bound = initial_radius ** 2

# Since K must be an integer, the maximum value is the integer part of the bound minus 1.
max_k = k_upper_bound - 1

# Output the explanation and the result.
print("The design requires that the slope of sightlines to the center increases with distance.")
print("This leads to the inequality: K < r^2")
print(f"The most restrictive case is for the smallest radius, r = {initial_radius} m.")
print(f"Substituting this value, we get the final equation: K < {initial_radius}^2")
print(f"This simplifies to: K < {k_upper_bound}")
print(f"Since K must be an integer, the maximum value it can take is {max_k}.")

print(f"\n<<<{{{max_k}}}>>>")