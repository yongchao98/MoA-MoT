# A great icosahedron is a non-convex polyhedron.
# Its outer-hull (or convex hull) is a regular icosahedron.
# A regular icosahedron is a Platonic solid with 20 faces, 12 vertices, and 30 edges.
# The question asks for the minimal number of triangles to triangulate this outer-hull.

# The faces of a regular icosahedron are already triangles.
# Therefore, the minimal triangulation is simply the set of its 20 faces.

num_triangles = 20

print("The outer-hull of the great icosahedron is a regular icosahedron.")
print("A regular icosahedron's surface is already composed of triangles.")
print("The minimal number of triangles for the triangulation is the number of its faces.")
print("The final number is:")
print(num_triangles)