import fractions

# Define the relationship between the outer and inner radii.
# Let's assume an inner radius of 1 for simplicity.
# This makes the outer radius 2, as per the problem statement (r_o = 2 * r_i).
r_i = 1
r_o = 2

# The formula for the ratio of maximum tangential stress to internal pressure is:
# Ratio = (r_o^2 + r_i^2) / (r_o^2 - r_i^2)

# Calculate the squared terms first
r_o_sq = r_o**2
r_i_sq = r_i**2

# Calculate the numerator and the denominator of the ratio formula
numerator = r_o_sq + r_i_sq
denominator = r_o_sq - r_i_sq

# Calculate the final ratio
ratio_decimal = numerator / denominator
ratio_fraction = fractions.Fraction(numerator, denominator)

# Print the calculation steps
print("The formula for the ratio of maximum tangential stress to internal pressure is:")
print("Ratio = (r_o^2 + r_i^2) / (r_o^2 - r_i^2)\n")
print(f"Given the outer radius is twice the inner radius, we can set r_i = {r_i} and r_o = {r_o}.")
print("Substituting these values into the formula:\n")
print(f"Ratio = ({r_o}^2 + {r_i}^2) / ({r_o}^2 - {r_i}^2)")
print(f"Ratio = ({r_o_sq} + {r_i_sq}) / ({r_o_sq} - {r_i_sq})")
print(f"Ratio = {numerator} / {denominator}\n")
print(f"The final ratio is {ratio_fraction}, which is approximately {ratio_decimal:.4f}.")