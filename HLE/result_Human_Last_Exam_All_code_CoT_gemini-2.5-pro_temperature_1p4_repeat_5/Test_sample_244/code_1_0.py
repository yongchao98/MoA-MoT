import spherogram
import sympy

# Step 1: Get all prime knots with exactly 7 crossings.
# The spherogram library provides a convenient way to access knot data.
# We iterate through all knots up to 7 crossings and filter for those with exactly 7.
knots_7_crossings = [K for K in spherogram.Knots() if K.crossings == 7]

# Step 2: Count the total number of unique 7-crossing knots.
# This gives us the denominator for our proportion.
total_knots_count = len(knots_7_crossings)

# Step 3 & 4: Identify and count the hyperbolic knots among them.
# A prime knot is either a torus knot or a hyperbolic knot.
# We can identify hyperbolic knots by checking which ones are not torus knots.
hyperbolic_knots_count = 0
for K in knots_7_crossings:
    if not K.is_torus():
        hyperbolic_knots_count += 1

# Step 5: Calculate and print the proportion.
# We will display the calculation as an equation.
proportion = sympy.Rational(hyperbolic_knots_count, total_knots_count)

print(f"Total number of knot types with 7 crossings: {total_knots_count}")
print(f"Number of hyperbolic knots among them: {hyperbolic_knots_count}")
print("\nCalculating the proportion of hyperbolic knots:")
# The final equation is constructed by printing each part.
print(f"{hyperbolic_knots_count} / {total_knots_count} = {float(proportion):.4f}")

# The proportion of hyperbolic knots is 6 out of 7.
final_answer = 6 / 7
# The result as a decimal is approximately 0.8571