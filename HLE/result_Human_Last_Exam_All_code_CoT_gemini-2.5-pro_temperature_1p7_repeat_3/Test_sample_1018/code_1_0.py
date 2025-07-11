import math

# Define the genera of the two surfaces
g1 = 31
g2 = 17

# Step 1: Compute the simplicial volume of each surface
# The formula for the simplicial volume of a surface Sigma_g with g >= 1 is ||Sigma_g|| = 4g - 4.
sv1 = 4 * g1 - 4
sv2 = 4 * g2 - 4

print(f"The simplicial volume of a surface Sigma_g is given by ||Sigma_g|| = 4g - 4 (for g >= 1).")
print(f"For Sigma_{g1}, the simplicial volume is 4 * {g1} - 4 = {sv1}")
print(f"For Sigma_{g2}, the simplicial volume is 4 * {g2} - 4 = {sv2}")
print("-" * 20)

# Step 2: Apply the product formula for simplicial volume
# The formula is ||M x N|| = C(dim(M)+dim(N), dim(M)) * ||M|| * ||N||
# Here, M = Sigma_31 and N = Sigma_17, both are 2-dimensional.
dim_m = 2
dim_n = 2
total_dim = dim_m + dim_n

# The binomial coefficient C(n, k)
coefficient = math.comb(total_dim, dim_m)

# The total simplicial volume
total_sv = coefficient * sv1 * sv2

print(f"The product formula is ||Sigma_{g1} x Sigma_{g2}|| = C({total_dim}, {dim_m}) * ||Sigma_{g1}|| * ||Sigma_{g2}||")
print(f"The binomial coefficient C({total_dim}, {dim_m}) is {coefficient}")
print("-" * 20)

# Step 3: Print the final calculation and result
print("Final Calculation:")
print(f"||\Sigma_{{{g1}}} x \Sigma_{{{g2}}}|| = {coefficient} * {sv1} * {sv2} = {total_sv}")
print("-" * 20)
print(f"The simplicial volume of Sigma_{g1} x Sigma_{g2} is {total_sv}.")
