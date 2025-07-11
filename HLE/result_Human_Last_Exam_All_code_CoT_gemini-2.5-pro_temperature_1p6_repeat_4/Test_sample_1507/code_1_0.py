# The outer-hull of the great icosahedron is a regular icosahedron.
# We need to find the number of its faces, which are already triangles.
# We can use Euler's formula for convex polyhedra: V - E + F = 2
# where V = vertices, E = edges, and F = faces.

# A regular icosahedron has 12 vertices.
V = 12
# A regular icosahedron has 30 edges.
E = 30
# The characteristic of any convex polyhedron is 2.
chi = 2

# Rearranging Euler's formula to solve for F: F = chi - V + E
F = chi - V + E

print(f"The calculation is based on Euler's formula for polyhedra: V - E + F = {chi}")
print(f"For the regular icosahedron (the outer-hull), V = {V} and E = {E}.")
print("The final equation to find the number of faces (F) is:")
print(f"{V} - {E} + F = {chi}")
print(f"Solving for F gives F = {chi} - {V} + {E}, which results in F = {F}.")
print(f"Therefore, the minimal number of triangles is {F}.")