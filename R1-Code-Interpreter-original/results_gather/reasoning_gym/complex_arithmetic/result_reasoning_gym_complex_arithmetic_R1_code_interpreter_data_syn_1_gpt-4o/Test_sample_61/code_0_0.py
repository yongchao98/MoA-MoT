# Given complex numbers
a, b = 61.0, -27.0
c, d = 1.0, -7.0

# Calculate the real part
real_part = (a * c + b * d) / (c**2 + d**2)

# Calculate the imaginary part
imaginary_part = (b * c - a * d) / (c**2 + d**2)

# Output the result
print(real_part, imaginary_part)