# Step 1: Define the total number of knot types with 7 crossings.
total_knots = 7

# Step 2: Identify the number of non-hyperbolic (torus) knots among them.
# The only non-hyperbolic prime knot with 7 crossings is the 7_1 knot, which is the (7,2)-torus knot.
non_hyperbolic_knots = 1

# Step 3: Calculate the number of hyperbolic knots.
hyperbolic_knots = total_knots - non_hyperbolic_knots

# Step 4: Calculate the proportion.
proportion = hyperbolic_knots / total_knots

# Step 5: Print the results, including the numbers used in the final equation.
print(f"Total number of knot types with 7 crossings: {total_knots}")
print(f"Number of hyperbolic knots with 7 crossings: {hyperbolic_knots}")
print(f"The proportion of hyperbolic knots is the ratio of hyperbolic knots to the total number of knots.")
print(f"Proportion = {hyperbolic_knots} / {total_knots} = {proportion}")
