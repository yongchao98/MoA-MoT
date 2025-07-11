import sympy

# Define symbols for the radii
# We don't need to define pressure as it will cancel out
a = sympy.Symbol('a')

# The outer radius 'b' is twice the inner radius 'a'
b = 2 * a

# The formula for the ratio of maximum tangential stress to internal pressure is:
# Ratio = (b^2 + a^2) / (b^2 - a^2)

# Calculate the numerator and denominator using the symbolic relationship
numerator = b**2 + a**2
denominator = b**2 - a**2

# Calculate the ratio
ratio = numerator / denominator

# Print the calculation steps
print(f"The formula for the ratio of maximum tangential stress to internal pressure is: (b^2 + a^2) / (b^2 - a^2)")
print(f"Given: Outer radius (b) = 2 * Inner radius (a)")
print(f"Substituting b = 2a into the formula:")
print(f"Ratio = ((2a)^2 + a^2) / ((2a)^2 - a^2)")
print(f"Ratio = (4a^2 + a^2) / (4a^2 - a^2)")

# Use sympy to evaluate and simplify the final expression
final_numerator = 4+1
final_denominator = 4-1
print(f"Ratio = ({final_numerator}a^2) / ({final_denominator}a^2)")
print(f"After cancelling a^2, the final ratio is {final_numerator}/{final_denominator}.")
print(f"As a decimal, this is approximately {float(ratio.evalf()):.3f}.")