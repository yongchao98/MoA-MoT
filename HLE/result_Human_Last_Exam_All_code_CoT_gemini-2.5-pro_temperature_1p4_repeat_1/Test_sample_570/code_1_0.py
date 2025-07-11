import math

# The minimal area is given by the formula 12 / (π^2).
# Here, π (pi) is the mathematical constant.

# Define the numbers in the equation
numerator = 12
pi_value = math.pi
exponent = 2

# Calculate the result
result = numerator / (pi_value ** exponent)

# As requested, output each number in the final equation, along with the result.
print(f"The final equation for the minimal area is: A = {numerator} / (π^{exponent})")
print(f"Numerator: {numerator}")
print(f"Value of π: {pi_value}")
print(f"Exponent: {exponent}")
print(f"Final Area A = {result}")