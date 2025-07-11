# Let's assign a nominal value to the inner radius 'a' to perform the calculation.
# Since we are calculating a ratio, the actual value of 'a' does not matter.
a = 1

# The outer radius 'b' is given as twice the inner radius 'a'.
b = 2 * a

# The formula for the ratio of the maximum tangential stress to the internal pressure is:
# Ratio = (b^2 + a^2) / (b^2 - a^2)

# Calculate the numerator and denominator of the equation using the given condition.
# The numerator corresponds to (b^2 + a^2)
# The denominator corresponds to (b^2 - a^2)
numerator_val = b**2 + a**2
denominator_val = b**2 - a**2

# The final numerical ratio
ratio_val = numerator_val / denominator_val

print("For a thick-walled cylinder with internal pressure, the ratio of maximum tangential stress to internal pressure is:")
print("Ratio = (b^2 + a^2) / (b^2 - a^2)")
print("\nGiven the condition that the outer radius (b) is twice the inner radius (a),")
print("we substitute b=2a, which simplifies the formula to 5a^2 / 3a^2.")
print("After canceling the common term (a^2), we get the final numerical equation:\n")
# Print the final equation with the calculated numbers
print(f"Ratio = {int(numerator_val)} / {int(denominator_val)}")
print(f"\nThis fraction, 5/3, is equivalent to approximately {ratio_val:.4f}.")