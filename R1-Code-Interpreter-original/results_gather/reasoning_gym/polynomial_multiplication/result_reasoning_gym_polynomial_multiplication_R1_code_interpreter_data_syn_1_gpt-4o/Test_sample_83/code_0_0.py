# Define the coefficients
a1, a2 = -37, -129
b1, b2 = 71, -69

# Calculate each term
term1 = a1 * b1  # Coefficient of y**6
term2 = a1 * b2  # Coefficient of y**3
term3 = a2 * b1  # Coefficient of y**5
term4 = a2 * b2  # Coefficient of y**2

# Print the results
print(f"{term1}*y**6 + {term3}*y**5 + {term2}*y**3 + {term4}*y**2")