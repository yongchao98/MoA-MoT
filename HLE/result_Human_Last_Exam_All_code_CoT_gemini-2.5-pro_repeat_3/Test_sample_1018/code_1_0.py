import math

def simplicial_volume_surface(g):
    """Calculates the simplicial volume of a surface of genus g."""
    if g < 1:
        return 0
    # For g >= 1, the formula is 4g - 4
    return 4 * g - 4

# Define the genera of the two surfaces
g1 = 31
g2 = 17

# Dimensions of the surfaces are both 2
dim1 = 2
dim2 = 2
total_dim = dim1 + dim2

# Step 1: Calculate the simplicial volume of the first surface
sv1 = simplicial_volume_surface(g1)

# Step 2: Calculate the simplicial volume of the second surface
sv2 = simplicial_volume_surface(g2)

# Step 3: Calculate the binomial coefficient C(4, 2)
binom_coeff = math.comb(total_dim, dim1)

# Step 4: Calculate the simplicial volume of the product
product_sv = binom_coeff * sv1 * sv2

# Print the final result showing the full equation
print(f"The simplicial volume of Sigma_{g1} x Sigma_{g2} is calculated as follows:")
print(f"||Sigma_{{{g1}}} x Sigma_{{{g2}}}|| = C({total_dim}, {dim1}) * ||Sigma_{{{g1}}}|| * ||Sigma_{{{g2}}}||")
print(f"Where ||Sigma_g|| = 4g - 4.")
print(f"||Sigma_{{{g1}}}|| = 4 * {g1} - 4 = {sv1}")
print(f"||Sigma_{{{g2}}}|| = 4 * {g2} - 4 = {sv2}")
print(f"C({total_dim}, {dim1}) = {binom_coeff}")
print("\nPutting it all together:")
print(f"||Sigma_{{{g1}}} x Sigma_{{{g2}}}|| = {binom_coeff} * {sv1} * {sv2} = {product_sv}")