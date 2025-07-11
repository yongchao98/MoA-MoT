import math

# Calculate the values of the Gamma functions
gamma_1_3 = math.gamma(1/3)
gamma_2_3 = math.gamma(2/3)

# Calculate the product
product_gamma = gamma_1_3 * gamma_2_3

# Calculate the simplified expression
simplified_value = 2 * math.pi / math.sqrt(3)

# Print the final equation with numerical values
print(f"The analytical value is derived from Gamma(1/3) * Gamma(2/3).")
print(f"Numerically:")
print(f"{gamma_1_3:.8f} * {gamma_2_3:.8f} = {product_gamma:.8f}")
print(f"This is equal to 2*pi/sqrt(3):")
print(f"2 * {math.pi:.8f} / {math.sqrt(3):.8f} = {simplified_value:.8f}")
print(f"Final analytical expression: Gamma(1/3) * Gamma(2/3) = 2*pi/sqrt(3)")