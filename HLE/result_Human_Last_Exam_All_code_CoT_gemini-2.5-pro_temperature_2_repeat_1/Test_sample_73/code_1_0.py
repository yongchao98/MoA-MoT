# Step 1: Define the Euler characteristic of the base region X.
# Based on the literature, the allowed configuration space on the parameter torus
# is topologically a disk with 4 holes.
num_holes = 4
# The Euler characteristic of a disk with k holes is 1 - k.
chi_X = 1 - num_holes
print(f"The base region X is a disk with {num_holes} holes.")
print(f"The Euler characteristic of X is 1 - {num_holes} = {chi_X}.")
print("-" * 20)

# Step 2: Define the Euler characteristic of the boundary of X.
# The boundary of a disk with 4 holes consists of 5 disjoint circles.
# The Euler characteristic of a circle is 0.
num_boundary_components = num_holes + 1
chi_boundary_X = 0
print(f"The boundary of X consists of {num_boundary_components} circles.")
print(f"The Euler characteristic of the boundary is {chi_boundary_X}.")
print("-" * 20)

# Step 3: Calculate the Euler characteristic of the configuration manifold M.
# The manifold M is a double cover of X branched over its boundary.
# The formula is chi(M) = 2 * chi(X) - chi(boundary_X).
chi_M = 2 * chi_X - chi_boundary_X
print("The configuration space M is a double cover of X.")
print(f"Its Euler characteristic is calculated as chi(M) = 2 * chi(X) - chi(boundary_X)")
print(f"chi(M) = 2 * ({chi_X}) - {chi_boundary_X} = {chi_M}")
print("-" * 20)

# Step 4: Calculate the genus g.
# For a closed, orientable surface, the formula is chi(M) = 2 - 2g.
# We solve for g: g = (2 - chi(M)) / 2.
genus = (2 - chi_M) / 2
print("The genus g is related to the Euler characteristic by chi(M) = 2 - 2g.")
print(f"Solving for g: g = (2 - chi(M)) / 2")
print(f"g = (2 - ({chi_M})) / 2 = {int(genus)}")
print("-" * 20)

print(f"The final genus of the surface is {int(genus)}.")