# Define the components of the complex numbers
a, b = 1.0, 8.0
c, d = -7.0, -10.0

# Calculate the real and imaginary parts
real_part = (a * c) - (b * d)
imaginary_part = (a * d) + (b * c)

# Output the result as a complex number
result = f"{real_part} + {imaginary_part}i"
print(result)