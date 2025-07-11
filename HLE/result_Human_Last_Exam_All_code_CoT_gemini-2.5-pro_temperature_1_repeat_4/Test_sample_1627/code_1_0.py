#
# This script calculates the ratio of maximum tangential stress to internal pressure
# for a thick-walled cylinder where the outer radius is twice the inner radius.
#

# Let's represent the inner radius 'a' and outer radius 'b' with numerical values
# that satisfy the condition b = 2a.
# We can set a = 1, which makes b = 2. The specific values don't matter, only their ratio.
a = 1
b = 2

# The formula for maximum tangential stress (sigma_t_max) in a thick-walled cylinder
# under internal pressure (P_i) is:
# sigma_t_max = P_i * (b^2 + a^2) / (b^2 - a^2)

# We want to find the ratio: sigma_t_max / P_i
# Ratio = (b^2 + a^2) / (b^2 - a^2)

# Let's calculate the numerator and denominator of this ratio.
numerator = b**2 + a**2
denominator = b**2 - a**2

# Now, let's print the equation with the calculated values.
print("The ratio of maximum tangential stress to internal pressure is given by the formula:")
print("Ratio = (b^2 + a^2) / (b^2 - a^2)")
print("\nSubstituting b = 2a (or using a=1, b=2):")
print(f"Ratio = ({b}^2 + {a}^2) / ({b}^2 - {a}^2)")
print(f"Ratio = ({b**2} + {a**2}) / ({b**2} - {a**2})")
print(f"Ratio = {numerator} / {denominator}")

# Calculate the final ratio
final_ratio_fraction = f"{int(numerator)}/{int(denominator)}"
final_ratio_decimal = numerator / denominator

print(f"\nThe final ratio is {final_ratio_fraction}, which is approximately {final_ratio_decimal:.4f}.")