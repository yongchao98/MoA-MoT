# Step 1: Define the components of the cellular decomposition of the surface.
# This decomposition is derived from viewing the configuration space as a 4-sheeted cover of a disk.

# F: The number of faces is the number of sheets in the cover.
F = 4

# E: The number of edges is determined by how the boundaries of the sheets are glued.
# The boundary of the base disk consists of two arcs. The four copies of each arc are identified in pairs.
E = 2 + 2

# V: The number of vertices is the number of points where the boundary arcs meet.
# At each of the two meeting points, all four sheets are joined.
V = 2

# Step 2: Calculate the Euler characteristic (chi) using the formula V - E + F.
chi = V - E + F

# Step 3: For a closed, orientable surface, the genus (g) is related to the
# Euler characteristic by the formula chi = 2 - 2g. We can solve for g.
# g = (2 - chi) / 2
g = (2 - chi) / 2

# Step 4: Print the full calculation and the final result.
print("To find the genus, we first calculate the Euler characteristic of the surface, chi = V - E + F.")
print(f"Number of Vertices (V): {V}")
print(f"Number of Edges (E): {E}")
print(f"Number of Faces (F): {F}")
print(f"Calculation: chi = {V} - {E} + {F} = {chi}")
print("\nNext, we find the genus (g) using the formula chi = 2 - 2g.")
print(f"Calculation: g = (2 - chi) / 2 = (2 - {chi}) / 2 = {int(g)}")
print(f"\nThe genus of the configuration space is {int(g)}.")