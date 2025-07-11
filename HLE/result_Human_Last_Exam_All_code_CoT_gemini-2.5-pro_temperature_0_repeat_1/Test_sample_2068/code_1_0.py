import math

# Given parameters
n = 2025
A = 10**15
B = 10**20

# The final expression to be calculated is A * B / (n - 1)
# as derived in the explanation.

# Numerator of the expression
numerator = A * B

# Denominator of the expression
denominator = n - 1

# Calculate the final result
result = numerator / denominator

# Output the numbers used in the final equation and the result
print(f"The value of A is: {A}")
print(f"The value of B is: {B}")
print(f"The value of n is: {n}")
print(f"The final equation is A * B / (n - 1)")
print(f"The value of the numerator (A * B) is: {numerator:.0e}")
print(f"The value of the denominator (n - 1) is: {denominator}")
print(f"The final result is: {result}")