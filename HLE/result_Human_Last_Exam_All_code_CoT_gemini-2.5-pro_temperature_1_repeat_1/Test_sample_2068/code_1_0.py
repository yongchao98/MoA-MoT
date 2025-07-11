import math

# Given constants
n = 2025
A = 10**15
B = 10**20

# The derivation shows that the final value is given by the formula (A * B) / (n - 1)
# Here we calculate the numerator and the denominator of the final expression
numerator = A * B
denominator = n - 1

# Calculate the final result
result = numerator / denominator

# Print the components of the calculation as requested
print("The value is calculated by the formula: (A * B) / (n - 1)")
print(f"A = {A}")
print(f"B = {B}")
print(f"n = {n}")
print(f"The equation for the value is: ({A} * {B}) / ({n} - 1)")
print(f"The numerator evaluates to: {numerator}")
print(f"The denominator evaluates to: {denominator}")
print(f"The final value is: {result}")