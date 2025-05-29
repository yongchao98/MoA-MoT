# Define the coefficients for each term
coefficients = {
    'z**5': -95 * 54,
    'z**4': -58 * 54,
    'z**3': -95 * 34,
    'z**2': -95 * (-14) - 58 * 34,
    'z': -58 * (-14)
}

# Calculate each term
z5 = coefficients['z**5']
z4 = coefficients['z**4']
z3 = coefficients['z**3']
z2 = coefficients['z**2']
z1 = coefficients['z']

# Print the results
print(f"{z5}*z**5 + {z4}*z**4 + {z3}*z**3 + {z2}*z**2 + {z1}*z")