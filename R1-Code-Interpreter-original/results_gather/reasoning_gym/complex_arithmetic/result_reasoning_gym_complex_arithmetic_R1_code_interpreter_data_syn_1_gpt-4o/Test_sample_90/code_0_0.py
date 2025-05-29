# Define the complex numbers
complex1 = (9.0, -1.0)  # (real, imaginary)
complex2 = (-8.0, 1.0)  # (real, imaginary)

# Subtract the real parts
real_part = complex1[0] - complex2[0]

# Subtract the imaginary parts
imaginary_part = complex1[1] - complex2[1]

# Print the result as a complex number
print(f"({real_part} + {imaginary_part}i)")