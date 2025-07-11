# Define the spatial dimensionality
d = 3

# For a scalar field theory, the number of components of the order parameter is n=1.
# The quantitative scaling of the specific heat exponent alpha (α) with
# dimensionality (d) is found using the epsilon expansion (ε = 4 - d).
# To first order, the relation for n=1 simplifies to: α ≈ (4 - d) / 6.

# Calculate the value of alpha for the given dimensionality
alpha_value = (4 - d) / 6

# We will now print the final equation, showing each number involved in the calculation.
print(f"The first-order scaling relation for α with dimensionality d is: α ≈ (4 - d) / 6")
print(f"For d = {d}:")
print(f"α ≈ (4 - {d}) / 6")
print(f"α ≈ {4-d} / 6")
print(f"α ≈ {alpha_value}")