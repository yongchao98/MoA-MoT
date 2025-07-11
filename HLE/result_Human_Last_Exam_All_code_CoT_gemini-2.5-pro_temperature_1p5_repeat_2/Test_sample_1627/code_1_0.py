# 1. Define variables based on the problem statement.
# The outer radius 'b' is twice the inner radius 'a'.
# We can use simple representative values since the ratio is dimensionless.
a = 1
b = 2

# 2. The formula for the ratio of maximum tangential stress to internal pressure is:
# Ratio = (a^2 + b^2) / (b^2 - a^2)

# 3. Calculate the components of the formula.
a_squared = a**2
b_squared = b**2
numerator = a_squared + b_squared
denominator = b_squared - a_squared

# 4. Calculate the final ratio.
ratio_value = numerator / denominator

# 5. Print the step-by-step calculation and the final equation.
print("This script calculates the ratio of maximum tangential stress to internal pressure for a thick-walled cylinder.")
print("The governing formula for the ratio is: (a^2 + b^2) / (b^2 - a^2)")
print(f"Given that the outer radius is twice the inner radius, we can set a = {a} and b = {b}.")
print("\nCalculating the equation step-by-step:")
print(f"Ratio = ({a}^2 + {b}^2) / ({b}^2 - {a}^2)")
print(f"Ratio = ({a_squared} + {b_squared}) / ({b_squared} - {a_squared})")
print(f"Ratio = {numerator} / {denominator}")
print(f"\nThe final ratio is 5/3, which is approximately {ratio_value:.4f}.")