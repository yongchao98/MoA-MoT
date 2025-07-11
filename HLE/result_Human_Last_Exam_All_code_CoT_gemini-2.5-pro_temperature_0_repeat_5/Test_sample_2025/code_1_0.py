# The value at which to evaluate the PDF
z = 0.2

# The formula for the PDF is f_Z(z) = 6 * z * (1 - z)
# We define the coefficients and terms for the final output
term1 = 6
term2 = z
term3 = 1 - z

# Calculate the result
result = term1 * term2 * term3

# Print the explanation, the equation with the numbers, and the final result.
print(f"The PDF of Z is given by the formula f_Z(z) = 6 * z * (1 - z).")
print(f"We want to calculate the value of this function at z = {z}.")
print(f"The calculation is:")
print(f"f({z}) = {term1} * {term2} * (1 - {term2})")
print(f"f({z}) = {term1} * {term2} * {term3}")
print(f"f({z}) = {result}")