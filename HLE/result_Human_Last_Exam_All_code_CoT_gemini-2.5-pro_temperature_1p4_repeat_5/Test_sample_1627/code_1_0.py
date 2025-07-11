import fractions

# Define the ratio of the outer radius to the inner radius
k = 2

# Calculate the numerator of the stress ratio formula: k^2 + 1
numerator = k**2 + 1

# Calculate the denominator of the stress ratio formula: k^2 - 1
denominator = k**2 - 1

# Calculate the final ratio
ratio_decimal = numerator / denominator
ratio_fraction = fractions.Fraction(numerator, denominator)

# Print the final equation with the calculated values
print(f"The ratio of maximum tangential stress to internal pressure is calculated as:")
print(f"Ratio = (k^2 + 1) / (k^2 - 1)")
print(f"Ratio = ({k}^2 + 1) / ({k}^2 - 1)")
print(f"Ratio = ({k**2} + 1) / ({k**2} - 1)")
print(f"Ratio = {int(numerator)} / {int(denominator)}")
print(f"\nThe ratio as a fraction is: {ratio_fraction}")
print(f"The ratio as a decimal is: {ratio_decimal:.4f}")