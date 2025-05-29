# Define the components of the complex numbers
a, b = -33.0, -21.0
c, d = -3.0, 3.0

# Calculate the real part of the result
real_part = (a * c + b * d)

# Calculate the imaginary part of the result
imaginary_part = (b * c - a * d)

# Calculate the denominator
denominator = (c ** 2 + d ** 2)

# Calculate the real and imaginary parts of the division
real_result = real_part / denominator
imaginary_result = imaginary_part / denominator

# Print the result as a complex number
print(f"{real_result} + {imaginary_result}i")