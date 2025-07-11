import math

# Given parameters
n = 2025
A = 10**15
B = 10**20

# The problem asks for the value of S / (2025 * pi).
# As derived in the plan, S = (n * pi * A * B) / (n - 1).
# So, the expression to calculate is ( (n * pi * A * B) / (n - 1) ) / (n * pi).
# The terms n and pi cancel out, leaving (A * B) / (n - 1).

# Calculate the numerator A * B
numerator = A * B

# Calculate the denominator n - 1
denominator = n - 1

# Calculate the final result
result = numerator / denominator

# Output the final equation with the given numbers
print(f"The final expression to calculate is (A * B) / (n - 1)")
print(f"Based on the given values, this becomes:")
final_equation_str = f"({A} * {B}) / ({n} - 1)"
print(final_equation_str)
print(f"= {numerator} / {denominator}")
print(f"Result: {result}")
