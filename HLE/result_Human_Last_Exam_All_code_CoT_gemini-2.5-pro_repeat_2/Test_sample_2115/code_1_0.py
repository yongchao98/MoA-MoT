import math

# Based on the analytical solution, we calculate the value of the integral.
# The formula is: -3 * ln((e^2 + e + 1) / 3)

# Define the constants from the derived formula
c1 = -3
c2 = 2
c3 = 1
c4 = 1
c5 = 3

# Calculate the term inside the logarithm
e_squared = math.exp(c2)
e_one = math.exp(c3)
numerator = e_squared + e_one + c4
denominator = c5
inner_expression = numerator / denominator

# Calculate the final result
result = c1 * math.log(inner_expression)

# Print the equation with the numbers to show the calculation steps
print(f"The value is calculated by the equation: {c1} * ln((e^{c2} + e^{c3} + {c4}) / {c5})")
print(f"The final result is: {result}")