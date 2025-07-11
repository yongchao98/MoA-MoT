#
#
# This script calculates the ratio of maximum tangential stress to internal pressure
# for a thick-walled cylinder where the outer radius is twice the inner radius.
#

# Let's define the radii for calculation. We can set the inner radius (ri) to 1
# for simplicity, which makes the outer radius (ro) equal to 2.
ri = 1
ro = 2

# The formula for the ratio is: (ro^2 + ri^2) / (ro^2 - ri^2)

# Calculate the numerator of the ratio formula
numerator = ro**2 + ri**2

# Calculate the denominator of the ratio formula
denominator = ro**2 - ri**2

# Calculate the final ratio
ratio = numerator / denominator

# Print the explanation and the final equation with numbers
print("The formula for the ratio of maximum tangential stress to internal pressure is:")
print("Ratio = (ro^2 + ri^2) / (ro^2 - ri^2)")
print(f"Given ro = {ro} * ri, we substitute the values (assuming ri = {ri}).")
print("\nThe final equation with the numbers substituted is:")
print(f"Ratio = ({ro}^2 + {ri}^2) / ({ro}^2 - {ri}^2)")
print(f"Ratio = ({ro**2} + {ri**2}) / ({ro**2} - {ri**2})")
print(f"Ratio = {numerator} / {denominator}")
print(f"\nThe calculated ratio is {numerator}/{denominator}, which is approximately {ratio:.4f}.")
