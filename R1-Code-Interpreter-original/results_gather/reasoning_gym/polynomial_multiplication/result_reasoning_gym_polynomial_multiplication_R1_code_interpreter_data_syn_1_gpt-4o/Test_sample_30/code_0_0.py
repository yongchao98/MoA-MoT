# Define the coefficients
coeff1 = 82
coeff2 = 21
coeff3 = -86
coeff4 = -31

# Calculate each term
term1 = coeff1 * coeff3  # Coefficient for z**5
term2 = coeff1 * coeff4  # Coefficient for z**3
term3 = coeff2 * coeff3  # Coefficient for z**3
term4 = coeff2 * coeff4  # Coefficient for z

# Combine like terms
z5_term = term1
z3_term = term2 + term3
z1_term = term4

# Print the simplified expression
print(f"Simplified expression: {z5_term}*z**5 + {z3_term}*z**3 + {z1_term}*z")