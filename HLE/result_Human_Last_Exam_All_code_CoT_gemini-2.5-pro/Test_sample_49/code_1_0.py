# Plan:
# 1. Define the f-vector for the base of the polytope, which is a 3D square pyramid.
# 2. Use the pyramid construction formulas to calculate the f-vector of the final 4-polytope.
# 3. Print the resulting f-vector and verify it with the Euler-Poincaré formula.

# Step 1: Define the f-vector of the base (a square pyramid).
# A square pyramid has 5 vertices, 8 edges, and 5 faces.
f0_base = 5
f1_base = 8
f2_base = 5

print(f"The base is a square pyramid with f-vector (vertices, edges, faces): ({f0_base}, {f1_base}, {f2_base})")
print("-" * 30)

# Step 2: Calculate the f-vector of the 4-polytope (pyramid over the square pyramid).
f0 = f0_base + 1          # Vertices
f1 = f1_base + f0_base      # Edges
f2 = f2_base + f1_base      # 2-Faces
f3 = f2_base + 1          # 3-Faces (Facets)

print("The f-vector of the resulting 4-polytope is (f0, f1, f2, f3).")
print(f"f0 (vertices) = {f0}")
print(f"f1 (edges)    = {f1}")
print(f"f2 (2-faces)  = {f2}")
print(f"f3 (facets)   = {f3}")
print("-" * 30)

# Step 3: Verify with the Euler-Poincaré formula for 4-polytopes (f0 - f1 + f2 - f3 = 0).
euler_char = f0 - f1 + f2 - f3
print("Checking the Euler-Poincaré formula for a 4-polytope:")
print(f"f0 - f1 + f2 - f3 = {f0} - {f1} + {f2} - {f3} = {euler_char}")

if euler_char == 0:
    print("The formula is satisfied.")
else:
    print("The formula is NOT satisfied.")

print("-" * 30)
print(f"The final f-vector is ({f0}, {f1}, {f2}, {f3}).")
