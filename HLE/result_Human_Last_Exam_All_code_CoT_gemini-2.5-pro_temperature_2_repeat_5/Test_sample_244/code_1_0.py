# Step 1: Identify the total number of knots with 7 crossings.
# The standard knot catalog (Alexander-Briggs notation) lists 7 distinct knots
# with 7 crossings: 7_1, 7_2, 7_3, 7_4, 7_5, 7_6, and 7_7.
total_knots = 7

# Step 2: Identify non-hyperbolic (torus) knots.
# According to Thurston's classification, a prime knot is either a torus knot or hyperbolic.
# All 7-crossing knots are prime. We check which of them are torus knots.
# The knot 7_1 is the torus knot T(7,2).
# The other six knots (7_2 to 7_7) are not torus knots.
non_hyperbolic_knots = 1

# Step 3: Calculate the number of hyperbolic knots.
# This is the total number of knots minus the non-hyperbolic ones.
hyperbolic_knots = total_knots - non_hyperbolic_knots

# Step 4: Calculate the proportion.
proportion = hyperbolic_knots / total_knots

# Print the final result, including the numbers used in the equation.
print(f"There are {total_knots} distinct knot types with 7 crossings.")
print(f"Among these, {non_hyperbolic_knots} is a torus knot, which is not hyperbolic.")
print(f"This means there are {hyperbolic_knots} hyperbolic knots.")
print(f"The proportion of hyperbolic knots is the number of hyperbolic knots divided by the total number of knots.")
print(f"Proportion = {hyperbolic_knots} / {total_knots}")
print(f"Result: {proportion}")
