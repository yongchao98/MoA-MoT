# The given parameters from the problem
n = 2025
A = 10**15
B = 10**20

# The final equation for the value is (A * B) / (n - 1)
# First, let's define the components of this equation
numerator = A * B
denominator = n - 1

# Calculate the result
result = numerator / denominator

# As requested, output the numbers in the final equation and the result
print(f"Numerator: {numerator}")
print(f"Denominator: {denominator}")
print(f"Result: {result}")