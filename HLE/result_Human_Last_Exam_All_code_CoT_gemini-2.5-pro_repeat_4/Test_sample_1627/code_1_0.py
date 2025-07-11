import fractions

# The problem asks for the ratio of the maximum tangential stress (σ_t_max)
# to the internal pressure (P_i) for a thick-walled cylinder.

# The formula for this ratio is:
# Ratio = (b^2 + a^2) / (b^2 - a^2)
# where 'a' is the inner radius and 'b' is the outer radius.

# We are given that the outer radius is twice the inner radius (b = 2a).
# For calculation, we can assume a=1, which makes b=2. The 'a' term
# cancels out in the ratio.
a = 1
b = 2

# Calculate the squares of the radii
a_squared = a**2
b_squared = b**2

# Calculate the numerator and denominator of the ratio formula
numerator = b_squared + a_squared
denominator = b_squared - a_squared

# Calculate the final ratio as a decimal
final_ratio_decimal = numerator / denominator

# To get the exact fraction
final_ratio_fraction = fractions.Fraction(numerator, denominator)


# Print the final equation with the numbers plugged in
print("The ratio is calculated using the formula: (b^2 + a^2) / (b^2 - a^2)")
print("With b=2 and a=1 (since b=2a), the equation becomes:")
print(f"Ratio = ({b_squared} + {a_squared}) / ({b_squared} - {a_squared})")
print(f"Ratio = {int(numerator)} / {int(denominator)}")
print(f"The value of the ratio is {final_ratio_fraction}, which is approximately {final_ratio_decimal:.3f}.")
