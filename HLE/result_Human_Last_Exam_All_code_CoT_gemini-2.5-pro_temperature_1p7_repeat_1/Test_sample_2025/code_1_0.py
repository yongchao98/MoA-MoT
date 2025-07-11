# The probability density function of Z is f_Z(z) = 6z^2 - 6z + 2
# We need to evaluate this function at z = 0.2.

z = 0.2

# Calculate the components of the polynomial
term1_coeff = 6
term1_val = z**2
term1 = term1_coeff * term1_val

term2_coeff = -6
term2_val = z
term2 = term2_coeff * term2_val

term3 = 2

# Calculate the final value of the PDF at z=0.2
f_z_at_0_2 = term1 + term2 + term3

# The problem asks us to output the equation with the numbers.
print(f"The PDF is f(z) = 6*z^2 - 6*z + 2")
print(f"Evaluating at z = 0.2:")
print(f"f(0.2) = {term1_coeff}*({z})^2 - {abs(term2_coeff)}*{z} + {term3}")
print(f"f(0.2) = {term1_coeff}*{term1_val} - {abs(term2_coeff)}*{term2_val} + {term3}")
print(f"f(0.2) = {term1} - {abs(term2)} + {term3}")
print(f"f(0.2) = {term1 + term2} + {term3}")
print(f"f(0.2) = {f_z_at_0_2}")