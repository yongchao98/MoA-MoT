import math

# The final expression simplifies to (3/2) * 10**(10/3) + 37/4.
# We will calculate this value and print the components of the equation.

# Define the numbers in the final equation
coeff1 = 3/2
base = 10
exponent = 10/3
term2 = 37/4

# Calculate the two main terms of the expression
term1_val = coeff1 * (base ** exponent)
term2_val = term2

# Calculate the final result
result = term1_val + term2_val

# Print the final equation with its components
print(f"The final simplified expression is: ({coeff1}) * {base}**({exponent}) + ({term2})")
print(f"The value of the first term is: {term1_val}")
print(f"The value of the second term is: {term2_val}")
print(f"The final result is: {result}")
