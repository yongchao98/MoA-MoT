# Define the coefficients of the polynomials
coeff1 = [-39, -93]
coeff2 = [61, 133, 43]

# Calculate the resulting coefficients after multiplication
x3_term = coeff1[0] * coeff2[0]
x2_term = coeff1[0] * coeff2[1] + coeff1[1] * coeff2[0]
x1_term = coeff1[0] * coeff2[2] + coeff1[1] * coeff2[1]
constant_term = coeff1[1] * coeff2[2]

# Print the simplified polynomial
print(f"{x3_term}*x**3 + {x2_term}*x**2 + {x1_term}*x + {constant_term}")