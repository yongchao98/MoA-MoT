# The paper is folded twice in each dimension, creating a 4x4 grid.
n = 4

# Step 1: Calculate the number of internal edges from the holes.
# The number of holes corresponds to the number of interior vertices in the grid.
# An n x n grid has an (n-1) x (n-1) grid of interior vertices.
num_internal_vertices = (n - 1) * (n - 1)
num_holes = num_internal_vertices

# Each hole is formed by cuts from the 4 squares meeting at an interior vertex,
# so each hole will have 4 edges.
edges_per_hole = 4
total_internal_edges = num_holes * edges_per_hole

print("Calculation for Internal Edges (Holes):")
print(f"Number of holes = (4 - 1) * (4 - 1) = {num_holes}")
print(f"Edges per hole = {edges_per_hole}")
print(f"Total internal edges = {num_holes} * {edges_per_hole} = {total_internal_edges}")
print("-" * 20)

# Step 2: Calculate the number of external edges (the outer boundary).
# The original square has 4 corners. Cutting them turns the square into an 8-sided shape.
base_edges_after_corner_cuts = 8

# The cuts also create nicks along the sides. The number of non-corner vertices on the
# boundary of an n x n grid is 4 * (n - 1).
num_side_vertices = 4 * (n - 1)
num_nicks = num_side_vertices

# Each nick adds 1 edge to the total count of the outer boundary.
added_edges_from_nicks = num_nicks
total_external_edges = base_edges_after_corner_cuts + added_edges_from_nicks

print("Calculation for External Edges (Boundary):")
print(f"Edges from cutting 4 corners of the square = {base_edges_after_corner_cuts}")
print(f"Number of nicks on sides = 4 * (4 - 1) = {num_nicks}")
print(f"Total external edges = {base_edges_after_corner_cuts} + {num_nicks} = {total_external_edges}")
print("-" * 20)

# Step 3: Calculate the total number of edges.
total_edges = total_internal_edges + total_external_edges

print("Total Number of Edges:")
print(f"Total Edges = Internal Edges + External Edges")
print(f"Total Edges = {total_internal_edges} + {total_external_edges} = {total_edges}")
print("-" * 20)
