# Define the real and imaginary parts of the complex numbers
a, b = 2.0, -9.0
c, d = 8.0, 8.0

# Calculate the real part of the product
real_part = a * c - b * d

# Calculate the imaginary part of the product
imaginary_part = a * d + b * c

# Print the result as a complex number
print(f"{real_part} + {imaginary_part}i")