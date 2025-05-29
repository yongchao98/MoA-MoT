# Define the complex numbers
complex1 = (8.0, 9.0)  # (real, imaginary)
complex2 = (-1.0, 4.0) # (real, imaginary)

# Subtract the real parts
real_part = complex1[0] - complex2[0]

# Subtract the imaginary parts
imaginary_part = complex1[1] - complex2[1]

# Form the resulting complex number
result = (real_part, imaginary_part)

# Print the result
print(f"Result: {result[0]} + {result[1]}i")