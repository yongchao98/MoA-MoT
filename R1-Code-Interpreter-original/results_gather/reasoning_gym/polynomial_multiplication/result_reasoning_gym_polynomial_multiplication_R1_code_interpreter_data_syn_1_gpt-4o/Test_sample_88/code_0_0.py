# Define the coefficients and powers
coefficients = [28, -82]
powers = [2, 1]
second_poly = [-9, -25]
second_powers = [3, 0]

# Calculate each term
term1 = coefficients[0] * second_poly[0] * z**(powers[0] + second_powers[0])
term2 = coefficients[0] * second_poly[1] * z**(powers[0] + second_powers[1])
term3 = coefficients[1] * second_poly[0] * z**(powers[1] + second_powers[0])
term4 = coefficients[1] * second_poly[1] * z**(powers[1] + second_powers[1])

# Combine the terms
simplified_expression = f"{term1} + {term3} + {term2} + {term4}"

# Print the simplified expression
print(simplified_expression)