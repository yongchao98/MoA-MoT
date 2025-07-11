# Step 1: Define the Euler characteristic of the boundary curve (∂M).
# The boundary curve has 2 self-intersections (nodes/vertices) and 4 edges connecting them.
V_boundary = 2
E_boundary = 4
chi_boundary = V_boundary - E_boundary
print(f"The boundary curve ∂M has {V_boundary} vertices and {E_boundary} edges.")
print(f"Euler characteristic of the boundary, χ(∂M) = {V_boundary} - {E_boundary} = {chi_boundary}")
print("-" * 20)

# Step 2: Define the Euler characteristic of the region M.
# From the topological analysis of how the boundary curve divides the torus:
# χ(M) + χ(M') = -4.
# We identified χ(M) as -3.
chi_M = -3
print(f"The Euler characteristic of the region M, χ(M) = {chi_M}")
print("-" * 20)

# Step 3: Calculate the Euler characteristic of the configuration space S.
# S is a double cover of M, with the two sheets joined at the boundary.
# The formula is χ(S) = 2 * χ(M) - χ(∂M).
chi_S = 2 * chi_M - chi_boundary
print("The Euler characteristic of the configuration space S is calculated as:")
print(f"χ(S) = 2 * χ(M) - χ(∂M)")
print(f"χ(S) = 2 * {chi_M} - ({chi_boundary}) = {chi_S}")
print("-" * 20)

# Step 4: Calculate the genus g from the Euler characteristic of S.
# The formula is χ(S) = 2 - 2g.
# Rearranging for g: g = (2 - χ(S)) / 2.
genus = (2 - chi_S) // 2
print("The genus g is calculated from χ(S) = 2 - 2g:")
print(f"{chi_S} = 2 - 2g")
print(f"2g = 2 - ({chi_S})")
print(f"2g = {2 - chi_S}")
print(f"g = {(2 - chi_S)} / 2 = {genus}")
print("-" * 20)

print(f"The genus of the configuration space is {genus}.")
