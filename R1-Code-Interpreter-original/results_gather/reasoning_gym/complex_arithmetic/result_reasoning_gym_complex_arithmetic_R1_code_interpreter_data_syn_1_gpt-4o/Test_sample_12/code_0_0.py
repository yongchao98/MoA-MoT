# Define the complex numbers
complex1 = (7.0, -7.0)  # (real, imaginary)
complex2 = (-5.0, 2.0)  # (real, imaginary)

# Subtract the real parts and the imaginary parts
real_part = complex1[0] - complex2[0]
imaginary_part = complex1[1] - complex2[1]

# Form the resulting complex number
result = (real_part, imaginary_part)

# Print the result
print(result)