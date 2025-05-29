# Define the complex numbers
a, b = 135.0, 20.0
c, d = 5.0, 10.0

# Calculate the real and imaginary parts of the result
real_part = (a * c + b * d) / (c**2 + d**2)
imaginary_part = (b * c - a * d) / (c**2 + d**2)

# Print the result as a complex number
print(f"{real_part} + {imaginary_part}i")