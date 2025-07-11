# Step 1: Define the f-vector of the base polytope (a triangular bipyramid).
# f0_base is the number of vertices (v).
# f1_base is the number of edges (e).
# f2_base is the number of 2-faces (f).
f0_base = 5  # Vertices
f1_base = 9  # Edges
f2_base = 6  # Faces

# Step 2: Calculate the f-vector of the 4-polytope (pyramid over the base).
# The 4-polytope is given to have 6 vertices.
f0 = f0_base + 1

# The number of edges (1-faces) is the sum of edges in the base
# plus the number of vertices in the base (as new edges to the apex).
f1 = f1_base + f0_base

# The number of 2-faces is the sum of 2-faces in the base
# plus the number of edges in the base (as new 2-faces with the apex).
f2 = f2_base + f1_base

# The number of 3-faces (facets) is the sum of 2-faces in the base
# (which form new facets with the apex) plus the base itself.
f3 = f2_base + 1

# Step 3: Print the final f-vector.
print("The polytope is the pyramid over the triangular bipyramid.")
print(f"Its f-vector (f0, f1, f2, f3) is calculated as follows:")
print(f"f0 = {f0_base} + 1 = {f0}")
print(f"f1 = {f1_base} + {f0_base} = {f1}")
print(f"f2 = {f2_base} + {f1_base} = {f2}")
print(f"f3 = {f2_base} + 1 = {f3}")
print("\nFinal f-vector:")
print(f"({f0}, {f1}, {f2}, {f3})")
