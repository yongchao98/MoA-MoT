import math

# The value of the sum is given by the formula: 1 - (pi^2 / 12) + (ln(2)^2 / 2)

# Calculate each term of the expression
term1 = 1.0
term2 = math.pi**2 / 12
term3 = (math.log(2))**2 / 2

# Calculate the final result
result = term1 - term2 + term3

# Print the equation with symbolic and numeric values
print("The final equation is derived as:")
print("S = 1 - (pi^2 / 12) + (ln(2)^2 / 2)")
print("\nCalculating the numerical values for each term:")
print(f"Term 1: {term1}")
print(f"Term 2 (pi^2 / 12): {term2}")
print(f"Term 3 (ln(2)^2 / 2): {term3}")

# Print the final equation with the computed numbers
print("\nSubstituting the numerical values into the equation:")
print(f"{term1} - {term2} + {term3} = {result}")

# Print the final answer
print(f"\nThe final value of the sum is: {result}")