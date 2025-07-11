import numpy as np

# Based on the analytical solution, the integral evaluates to:
# -3 * (ln(e^2 + e + 1) - ln(3))
# Let's define the components of this final equation.

# The coefficient
coeff = -3

# The first term inside the logarithm
term1 = np.exp(2) + np.exp(1) + 1

# The second term inside the logarithm
term2 = 3

# Calculate the final result
result = coeff * (np.log(term1) - np.log(term2))

# Print the components of the final equation as requested
print("The definite integral evaluates to the expression: A * (ln(B) - ln(C))")
print(f"where the numbers in the equation are:")
print(f"A = {coeff}")
print(f"B = e^2 + e + 1 = {term1}")
print(f"C = {term2}")
print("\n---")

# Print the final numerical answer
print(f"The final calculated value of the spatial average is: {result}")