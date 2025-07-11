# This script calculates the ratio of maximum tangential stress to internal pressure
# for a thick-walled cylinder where the outer radius is twice the inner radius.

# Let 'k' be the ratio of the outer radius to the inner radius (ro/ri).
k = 2

# The formula for the ratio of maximum tangential stress to internal pressure is:
# Ratio = (ro^2 + ri^2) / (ro^2 - ri^2)
# By dividing the numerator and denominator by ri^2, this simplifies to:
# Ratio = (k^2 + 1) / (k^2 - 1)

# Calculate the intermediate values for demonstrating the equation
k_squared = k**2
numerator = k_squared + 1
denominator = k_squared - 1

# Print the calculation steps
print("The formula for the ratio of maximum tangential stress to internal pressure is (k^2 + 1) / (k^2 - 1), where k = outer_radius / inner_radius.")
print(f"For this problem, k = {k}.")
print("\nCalculation steps:")
print(f"Ratio = ({k}^2 + 1) / ({k}^2 - 1)")
print(f"      = ({k_squared} + 1) / ({k_squared} - 1)")
print(f"      = {int(numerator)} / {int(denominator)}")

print(f"\nThus, the final ratio of the maximum tangential stress to the internal pressure is {int(numerator)}/{int(denominator)}.")
