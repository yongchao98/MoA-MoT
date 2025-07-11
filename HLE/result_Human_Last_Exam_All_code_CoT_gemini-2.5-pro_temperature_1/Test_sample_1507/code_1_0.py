# The problem asks for the number of triangles on the outer-hull of the great icosahedron.
# Step 1: The outer-hull (or convex hull) of a great icosahedron is a regular icosahedron.
# Step 2: The number of visible triangles is the number of faces of this regular icosahedron.
# Step 3: We can use Euler's polyhedron formula (V - E + F = 2) to find the number of faces (F).

# Step 4: A regular icosahedron has 12 vertices (V) and 30 edges (E).
V = 12
E = 30

# Step 5: Calculate the number of faces (F) using the formula F = 2 - V + E.
F = 2 - V + E

print("The convex hull of the great icosahedron is a regular icosahedron.")
print("A regular icosahedron has V = {} vertices and E = {} edges.".format(V, E))
print("Using Euler's formula V - E + F = 2, we can find the number of faces (F).")
print("F = 2 - V + E")
print("F = 2 - {} + {}".format(V, E))
print("The minimal number of triangles needed is: {}".format(F))