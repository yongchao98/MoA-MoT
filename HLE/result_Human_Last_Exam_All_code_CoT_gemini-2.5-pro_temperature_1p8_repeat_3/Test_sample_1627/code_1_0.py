# Define the ratio of the outer radius to the inner radius
# The problem states r_o = 2 * r_i, so k = r_o / r_i = 2
k = 2.0

# The formula for the ratio of maximum tangential stress (sigma_t_max) to internal pressure (P_i)
# for a thick-walled cylinder is: (r_o^2 + r_i^2) / (r_o^2 - r_i^2)
# This can be simplified by dividing the numerator and denominator by r_i^2 to:
# (k^2 + 1) / (k^2 - 1)

# Calculate the numerator and denominator of the simplified formula
numerator = k**2 + 1
denominator = k**2 - 1

# Calculate the final ratio as a floating-point number
ratio_value = numerator / denominator

print("The formula for the ratio of maximum tangential stress to internal pressure is:")
print("Ratio = (r_o^2 + r_i^2) / (r_o^2 - r_i^2)")
print("\nGiven r_o = 2 * r_i, we can substitute this into the equation:")
print("Ratio = ((2*r_i)^2 + r_i^2) / ((2*r_i)^2 - r_i^2)")
print("Ratio = (4*r_i^2 + r_i^2) / (4*r_i^2 - r_i^2)")
print("Ratio = (5*r_i^2) / (3*r_i^2)")
print("\nThe r_i^2 terms cancel, which gives us the final equation for the ratio:")

# Print each number in the final equation
print(f"Ratio = {int(numerator)} / {int(denominator)}")

# Also print the decimal representation for clarity
print(f"\nThe numerical value is {ratio_value}, which is approximately {ratio_value:.3f}.")
