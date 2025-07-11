# 1. Define the total number of knot types with 7 crossings.
# As per the Rolfsen knot table, there are 7 such knots.
total_knots = 7

# 2. Define the number of non-hyperbolic knots with 7 crossings.
# Among the 7-crossing knots, only 7_1 is a torus knot, which is not hyperbolic.
non_hyperbolic_knots = 1

# 3. Calculate the number of hyperbolic knots.
# This is the total number of knots minus the non-hyperbolic ones.
hyperbolic_knots = total_knots - non_hyperbolic_knots

# 4. Calculate the proportion of hyperbolic knots.
proportion = hyperbolic_knots / total_knots

# As requested, print the numbers that form the final equation for the proportion.
print(f"Total number of 7-crossing knots: {total_knots}")
print(f"Number of hyperbolic 7-crossing knots: {hyperbolic_knots}")
print(f"The proportion of hyperbolic knots is calculated as follows:")
print(f"{hyperbolic_knots} / {total_knots}")

# For clarity, also print the resulting decimal value.
print(f"Result: {proportion}")