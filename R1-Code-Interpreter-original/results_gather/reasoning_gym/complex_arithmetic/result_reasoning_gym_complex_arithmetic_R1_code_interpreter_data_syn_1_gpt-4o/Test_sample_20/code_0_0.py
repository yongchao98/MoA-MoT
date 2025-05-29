# Define the complex numbers
numerator = complex(-51.0, 21.0)
denominator = complex(5.0, 1.0)

# Calculate the conjugate of the denominator
conjugate_denominator = complex(denominator.real, -denominator.imag)

# Multiply the numerator by the conjugate of the denominator
numerator_conjugate_product = numerator * conjugate_denominator

# Multiply the denominator by its conjugate
denominator_conjugate_product = denominator * conjugate_denominator

# Perform the division
result = numerator_conjugate_product / denominator_conjugate_product

# Print the result
print(result)