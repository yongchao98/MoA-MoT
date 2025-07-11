# Define the relationship between the radii.
# We can use arbitrary values that satisfy the condition, e.g., r_i = 1, r_o = 2.
ri = 1
ro = 2

# The formula for the ratio of maximum tangential stress to internal pressure is:
# Ratio = (ro^2 + ri^2) / (ro^2 - ri^2)

# Calculate the numerator and denominator
numerator = ro**2 + ri**2
denominator = ro**2 - ri**2

# Calculate the final ratio
ratio_decimal = numerator / denominator
ratio_fraction_str = f"{int(numerator)}/{int(denominator)}"

# Print the equation with the calculated values
print(f"The ratio of maximum tangential stress to internal pressure is calculated as follows:")
print(f"Ratio = (ro^2 + ri^2) / (ro^2 - ri^2)")
print(f"Ratio = ({ro**2} + {ri**2}) / ({ro**2} - {ri**2})")
print(f"Ratio = {int(numerator)} / {int(denominator)}")
print(f"The final ratio is {ratio_fraction_str}, which is approximately {ratio_decimal:.3f}.")
