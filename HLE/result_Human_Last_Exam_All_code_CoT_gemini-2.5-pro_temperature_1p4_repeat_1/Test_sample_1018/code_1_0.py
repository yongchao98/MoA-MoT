import math

# Step 1: Define the genera of the two surfaces.
g1 = 31
g2 = 17

print(f"The task is to compute the simplicial volume of the product of two surfaces: Sigma_{g1} and Sigma_{g2}.")

# Step 2: Calculate the simplicial volume for each surface.
# The formula for a surface with genus g >= 2 is ||Sigma_g|| = 4g - 4.
sv1 = 4 * g1 - 4
print(f"The simplicial volume of the first surface, ||Sigma_{g1}||, is 4 * {g1} - 4 = {sv1}.")

sv2 = 4 * g2 - 4
print(f"The simplicial volume of the second surface, ||Sigma_{g2}||, is 4 * {g2} - 4 = {sv2}.")

# Step 3: Apply the product formula.
# The formula is ||M x N|| = C(dim(M)+dim(N), dim(M)) * ||M|| * ||N||.
# The dimension of a surface is 2.
dim1 = 2
dim2 = 2
total_dim = dim1 + dim2

# Calculate the binomial coefficient C(4, 2).
coefficient = math.comb(total_dim, dim1)
print(f"The product formula requires the binomial coefficient C({dim1}+{dim2}, {dim1}) = C({total_dim}, {dim1}) = {coefficient}.")

# Step 4: Compute the final result.
final_volume = coefficient * sv1 * sv2
print(f"The simplicial volume of the product manifold Sigma_{g1} x Sigma_{g2} is:")
print(f"{coefficient} * {sv1} * {sv2} = {final_volume}")