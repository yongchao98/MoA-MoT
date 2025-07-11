import math

# The final simplified expression to calculate is (3/2) * 10^(10/3) + 37/4

# Define the components of the final equation
term1_coeff = 3/2
term1_base = 10
term1_exponent = 10/3
term2 = 37/4

# Calculate the value of the first term
term1_value = term1_coeff * (term1_base ** term1_exponent)

# Calculate the final result
result = term1_value + term2

# Print the components of the final equation as requested
print("The final calculation is based on the equation: (3/2) * 10^(10/3) + 37/4")
print(f"Value of the first coefficient (3/2): {term1_coeff}")
print(f"Value of the base: {term1_base}")
print(f"Value of the exponent (10/3): {term1_exponent}")
print(f"Value of the second term (37/4): {term2}")
print(f"Final result: {result}")
