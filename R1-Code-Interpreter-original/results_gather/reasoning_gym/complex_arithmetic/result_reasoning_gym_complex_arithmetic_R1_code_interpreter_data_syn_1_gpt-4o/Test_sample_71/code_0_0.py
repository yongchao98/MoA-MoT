# Define the complex numbers
numerator = complex(4.0, 68.0)
denominator = complex(-8.0, 9.0)

# Calculate the conjugate of the denominator
conjugate_denominator = complex(denominator.real, -denominator.imag)

# Perform the division
result = numerator * conjugate_denominator / (denominator * conjugate_denominator)

# Output the result
print(result)