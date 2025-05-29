# Define the complex numbers
complex1 = (10.0, -3.0)  # (real, imaginary)
complex2 = (-2.0, 5.0)   # (real, imaginary)

# Subtract the real parts and the imaginary parts
real_part = complex1[0] - complex2[0]
imaginary_part = complex1[1] - complex2[1]

# Form the resulting complex number
result = (real_part, imaginary_part)

# Print the result
print(f"Result: {result[0]} {'+' if result[1] >= 0 else '-'} {abs(result[1])}i")