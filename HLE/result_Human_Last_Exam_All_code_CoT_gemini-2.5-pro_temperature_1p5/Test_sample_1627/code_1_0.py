import fractions

# Problem statement gives the relationship between the outer and inner radius.
# Let's set the inner radius (r_i) to 1 for simplicity.
# This means the outer radius (r_o) is 2.
r_i = 1
r_o = 2

# The formula for the ratio of maximum tangential stress to internal pressure is:
# Ratio = (r_o^2 + r_i^2) / (r_o^2 - r_i^2)

# Calculate the squares of the radii
ro_sq = r_o**2
ri_sq = r_i**2

# Calculate the numerator and denominator of the ratio formula
numerator = ro_sq + ri_sq
denominator = ro_sq - ri_sq

# The final code should output each number in the final equation.
print("The ratio is calculated using the formula: (r_o^2 + r_i^2) / (r_o^2 - r_i^2)")
print("With r_o = 2 and r_i = 1, the equation becomes:")
print(f"({ro_sq} + {ri_sq}) / ({ro_sq} - {ri_sq})")
print(f"= {numerator} / {denominator}")

# Represent the result as a fraction
final_fraction = fractions.Fraction(numerator, denominator)
print(f"\nThe final ratio is {final_fraction}, which is approximately {float(final_fraction):.3f}.")