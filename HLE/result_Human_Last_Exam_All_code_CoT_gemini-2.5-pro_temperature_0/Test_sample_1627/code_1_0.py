import fractions

# Step 1 & 2: Define the formula for the ratio.
# The ratio of maximum tangential stress to internal pressure for a thick-walled cylinder is:
# Ratio = (b^2 + a^2) / (b^2 - a^2)
# where 'a' is the inner radius and 'b' is the outer radius.

# Step 3: Apply the given condition.
# The problem states that the outer radius is twice the inner radius (b = 2a).
# For calculation purposes, we can assume a = 1, which makes b = 2.
a = 1
b = 2

# Step 4: Calculate the terms of the equation.
a_squared = a**2
b_squared = b**2

numerator = b_squared + a_squared
denominator = b_squared - a_squared

# Calculate the final ratio
ratio_decimal = numerator / denominator
# For a more precise representation, calculate as a fraction
ratio_fraction = fractions.Fraction(numerator, denominator)


# Print the results, showing the calculation step-by-step.
print("The formula for the ratio is: (b^2 + a^2) / (b^2 - a^2)")
print(f"Given b = 2a, we can substitute a={a} and b={b}.")
print(f"The equation becomes: ({b}^2 + {a}^2) / ({b}^2 - {a}^2)")
print(f"= ({b_squared} + {a_squared}) / ({b_squared} - {a_squared})")
print(f"= {numerator} / {denominator}")
print(f"The resulting ratio is {ratio_decimal}, which is equal to the fraction {ratio_fraction}.")
