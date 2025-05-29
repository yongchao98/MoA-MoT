# Given complex numbers
a, b = -56.0, 34.0
c, d = 7.0, 3.0

# Calculate the real part
real_part = (a * c + b * d)

# Calculate the imaginary part
imaginary_part = (b * c - a * d)

# Calculate the denominator
denominator = (c ** 2 + d ** 2)

# Divide the real and imaginary parts by the denominator
real_result = real_part / denominator
imaginary_result = imaginary_part / denominator

# Print the result
print(f"({real_result}) + ({imaginary_result})i")