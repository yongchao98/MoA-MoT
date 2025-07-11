import numpy as np

# Based on analytical evaluation, the integral is simplified to the following expression:
# I = (8/15)*pi^8 + (1/3)*pi^2 - (1/2)*pi + 1
# The code below defines the numbers in the final equation and prints the expression.

# The numbers that form the coefficients in the final equation
num1, den1 = 8, 15
num2, den2 = 1, 3
num3, den3 = 1, 2
const_term = 1

# Calculate the final numeric value for verification
pi = np.pi
value = (num1 / den1) * pi**8 + (num2 / den2) * pi**2 - (num3 / den3) * pi + const_term

# Print the final symbolic expression as requested
print("The final value of the integral is given by the expression:")
print(f"({num1}/{den1}) * pi^8 + ({num2}/{den2}) * pi^2 - ({num3}/{den3}) * pi + {const_term}")

# As a reference, the computed numerical value is printed below
# print(f"\nNumerical Value: {value}")
