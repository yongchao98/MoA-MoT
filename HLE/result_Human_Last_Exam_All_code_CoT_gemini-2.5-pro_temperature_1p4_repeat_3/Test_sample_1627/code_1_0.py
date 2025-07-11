import fractions

# Step 1: Define the relationship between the radii.
# Let's set the inner radius (ri) to 1 for simplicity.
# The outer radius (ro) is twice the inner radius.
ri = 1
ro = 2 * ri

# Step 2: Use the formula for the ratio of maximum tangential stress to internal pressure.
# Ratio = (ro^2 + ri^2) / (ro^2 - ri^2)

# Calculate the squares of the radii
ro_sq = ro**2
ri_sq = ri**2

# Calculate the numerator and denominator of the ratio formula
numerator = ro_sq + ri_sq
denominator = ro_sq - ri_sq

# Step 3: Print the equation with the calculated values
# The prompt requires printing each number in the final equation.
print(f"The ratio σ_t_max / P_i is given by the formula: (ro² + ri²) / (ro² - ri²)")
print(f"Substituting ro = {ro} and ri = {ri}:")
print(f"Ratio = ({ro}² + {ri}²) / ({ro}² - {ri}²)")
print(f"      = ({ro_sq} + {ri_sq}) / ({ro_sq} - {ri_sq})")
print(f"      = {numerator} / {denominator}")

# Step 4: Calculate the final ratio as a decimal and a fraction
final_ratio_decimal = numerator / denominator
final_ratio_fraction = fractions.Fraction(numerator, denominator)

print(f"\nThe final value of the ratio is {final_ratio_decimal:.4f}, or {final_ratio_fraction.numerator}/{final_ratio_fraction.denominator} as a fraction.")