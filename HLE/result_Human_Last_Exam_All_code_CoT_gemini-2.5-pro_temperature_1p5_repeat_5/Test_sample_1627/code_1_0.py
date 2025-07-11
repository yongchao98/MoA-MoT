# Define the radii based on the problem statement.
# Let's set the inner radius 'a' to 1 for simplicity.
# The outer radius 'b' is twice the inner radius.
a = 1
b = 2 * a

# The formula for the ratio of maximum tangential stress (sigma_t_max)
# to internal pressure (P_i) in a thick-walled cylinder is:
# Ratio = (b^2 + a^2) / (b^2 - a^2)

# Calculate the numerator and the denominator of the formula
numerator = b**2 + a**2
denominator = b**2 - a**2

# Calculate the final ratio
ratio = numerator / denominator

# Print the calculation steps
print("The formula for the ratio of maximum tangential stress to internal pressure is:")
print("Ratio = (b^2 + a^2) / (b^2 - a^2)")
print(f"Given b = 2 * a, let a = {a}, then b = {b}.")
print("\nSubstituting these values into the formula:")
print(f"Ratio = ({b}^2 + {a}^2) / ({b}^2 - {a}^2)")
print(f"Ratio = ({b**2} + {a**2}) / ({b**2} - {a**2})")
print(f"Ratio = {numerator} / {denominator}")
print(f"\nThe calculated ratio is {ratio:.4f}, which is equivalent to the fraction {numerator}/{denominator}.")