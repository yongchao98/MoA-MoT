import math

# Based on the derivation, the final expression to evaluate is:
# (3/2) * 10^(10/3) + 37/4

# Define the components of the final equation
numerator1 = 3
denominator1 = 2
base = 10
exponent_num = 10
exponent_den = 3
numerator2 = 37
denominator2 = 4

# Calculate the terms
term1 = (numerator1 / denominator1) * (base ** (exponent_num / exponent_den))
term2 = numerator2 / denominator2

# Calculate the final result
result = term1 + term2

# Output the equation and the final answer
print("The final expression simplifies to the equation:")
print(f"({numerator1}/{denominator1}) * {base}^({exponent_num}/{exponent_den}) + ({numerator2}/{denominator2})")
print("\nCalculating the value:")
print(f"First term: {term1}")
print(f"Second term: {term2}")
print(f"Result = {term1} + {term2} = {result}")
