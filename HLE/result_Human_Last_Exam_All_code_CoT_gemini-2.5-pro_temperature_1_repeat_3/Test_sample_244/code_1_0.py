# Step 1: Define the list of all prime knots with 7 crossings.
# There are 7 such knots according to the Rolfsen knot table.
knots_7_crossings = ['7_1', '7_2', '7_3', '7_4', '7_5', '7_6', '7_7']
total_knots_count = len(knots_7_crossings)

# Step 2: Identify the non-hyperbolic knots among them.
# A prime knot is non-hyperbolic if it is a torus knot.
# For 7-crossing knots, only 7_1 (the torus knot T(2,7)) is non-hyperbolic.
non_hyperbolic_knots = ['7_1']
non_hyperbolic_knots_count = len(non_hyperbolic_knots)

# Step 3: Calculate the number of hyperbolic knots.
hyperbolic_knots_count = total_knots_count - non_hyperbolic_knots_count

# Step 4: Calculate the proportion of hyperbolic knots.
proportion = hyperbolic_knots_count / total_knots_count

# Print the final equation and the result.
print(f"Total number of 7-crossing knots: {total_knots_count}")
print(f"Number of hyperbolic 7-crossing knots: {hyperbolic_knots_count}")
print(f"Proportion of hyperbolic knots = {hyperbolic_knots_count} / {total_knots_count} = {proportion}")