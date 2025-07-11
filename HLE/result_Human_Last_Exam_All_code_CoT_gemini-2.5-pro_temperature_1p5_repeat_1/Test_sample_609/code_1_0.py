import math

# The problem asks for a general solution, which we derived as a formula.
# This script calculates the ratio for the example case where n=3.
# This corresponds to constructing a triangle (n=3) from a hexagon (2n=6).

# Value of n for the case of an equilateral triangle
n = 3

print(f"Calculating the area ratio for a {n}-sided polygon constructed from a {2*n}-sided polygon.")
print(f"The general formula for the ratio is: (1 + sec(\u03C0/n)) / 2")
print("-" * 20)

# Perform the calculation step-by-step
pi_val = math.pi
denominator = 2

# Calculate the argument for the trigonometric functions
# This is the first number in the equation: n
secant_argument = pi_val / n
print(f"Step 1: Calculate the argument for the secant function, \u03C0/n = {pi_val:.4f} / {n} = {secant_argument:.4f} radians.")

# Calculate the cosine
cos_val = math.cos(secant_argument)
print(f"Step 2: Calculate the cosine of this value, cos({secant_argument:.4f}) = {cos_val:.4f}.")

# Calculate the secant
secant_val = 1 / cos_val
print(f"Step 3: Calculate the secant (1/cosine), sec({secant_argument:.4f}) = {secant_val:.4f}.")

# Calculate the numerator of the formula
# These are the numbers 1 and the result from Step 3
numerator = 1 + secant_val
print(f"Step 4: Calculate the numerator of the formula, 1 + {secant_val:.4f} = {numerator:.4f}.")

# Calculate the final ratio
# These are the numbers from Step 4 and the denominator, 2
ratio = numerator / denominator
print(f"Step 5: Calculate the final ratio, {numerator:.4f} / {denominator} = {ratio:.1f}.")

print("-" * 20)
print(f"Final Answer: The area of the triangle is {ratio} times larger than the hexagon.")
