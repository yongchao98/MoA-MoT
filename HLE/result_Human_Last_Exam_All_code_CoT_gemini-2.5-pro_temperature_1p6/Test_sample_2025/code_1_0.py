# The value z for which we want to calculate the pdf
z = 0.2

# The PDF f_Z(z) for z in [0, 0.5] is given by the formula 3 - 12*z + 12*z^2.
# Since 0.2 is in this range, we can directly apply the formula.

# Define the coefficients of the polynomial
c0 = 3
c1 = -12
c2 = 12

# Calculate the value of the pdf at z = 0.2
fz_at_0_2 = c2 * z**2 + c1 * z + c0

# We print the calculation step by step
print(f"The PDF is calculated using the formula: f(z) = {c2}*z^2 + {c1}*z + {c0}")
print(f"For z = {z}:")
term1 = c2 * z**2
term2 = c1 * z
term3 = c0
print(f"f({z}) = {c2}*({z}^2) + ({c1})*({z}) + {c0}")
print(f"f({z}) = {c2}*{z**2} + ({term2}) + {c0}")
print(f"f({z}) = {term1} + ({term2}) + {c0}")
result = term1 + term2 + term3
print(f"f({z}) = {result}")

# The exact value can be represented as a fraction
# 1.08 = 108/100 = 27/25
print(f"The exact value is 27/25.")

# Final check of the calculation
# 3 - 12*(0.2) + 12*(0.04) = 3 - 2.4 + 0.48 = 0.6 + 0.48 = 1.08
# As a fraction: 3 - 12/5 + 12/25 = 75/25 - 60/25 + 12/25 = 27/25
# 27/25 = 1.08
print(f"Final calculated value is {result}, which is equal to 27/25.")