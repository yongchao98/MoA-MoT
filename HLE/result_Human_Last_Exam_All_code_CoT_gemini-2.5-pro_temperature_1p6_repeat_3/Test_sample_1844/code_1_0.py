# Plan:
# 1. Define the grid size based on the number of folds.
#    - Two pairs of folds (top-to-bottom, left-to-right) create a 4x4 grid.
# 2. Calculate the number of vertices of each type: corner, edge, and interior.
# 3. For each vertex type, determine how many edges are created by the cuts.
# 4. Sum the edges from all vertex types to get the total number of edges.

# n is the number of subdivisions along one side of the paper.
# Two horizontal and two vertical folds lead to n=4.
n = 4

# --- Calculations for each vertex type ---

# 1. Corner vertices
# There are always 4 corners on a square.
num_corner_vertices = 4
# Each corner cut creates 1 new edge on the final shape's boundary.
edges_per_corner = 1
total_edges_from_corners = num_corner_vertices * edges_per_corner

print(f"Number of corner vertices: {num_corner_vertices}")
print(f"Edges created per corner cut: {edges_per_corner}")
print(f"Subtotal of edges from corner cuts: {num_corner_vertices} * {edges_per_corner} = {total_edges_from_corners}\n")

# 2. Edge vertices (not including corners)
# Each of the 4 sides has (n-1) vertices between the corners.
num_edge_vertices = 4 * (n - 1)
# At an edge vertex, two small squares meet. Their cuts merge to form a V-shaped notch with 2 edges.
edges_per_edge_notch = 2
total_edges_from_edges = num_edge_vertices * edges_per_edge_notch

print(f"Number of edge vertices: 4 * ({n} - 1) = {num_edge_vertices}")
print(f"Edges created per edge notch: {edges_per_edge_notch}")
print(f"Subtotal of edges from edge notches: {num_edge_vertices} * {edges_per_edge_notch} = {total_edges_from_edges}\n")

# 3. Interior vertices
# The interior vertices form a grid of size (n-1)x(n-1).
num_interior_vertices = (n - 1) * (n - 1)
# At an interior vertex, four small squares meet. Their cuts merge to form a hole with 4 edges (a diamond/square shape).
edges_per_interior_hole = 4
total_edges_from_interior = num_interior_vertices * edges_per_interior_hole

print(f"Number of interior vertices: ({n} - 1) * ({n} - 1) = {num_interior_vertices}")
print(f"Edges created per interior hole: {edges_per_interior_hole}")
print(f"Subtotal of edges from interior holes: {num_interior_vertices} * {edges_per_interior_hole} = {total_edges_from_interior}\n")

# --- Final Calculation ---

# The total number of edges is the sum of the edges from all three types of locations.
total_edges = total_edges_from_corners + total_edges_from_edges + total_edges_from_interior

print(f"Total number of edges = (edges from corners) + (edges from edge notches) + (edges from interior holes)")
print(f"Total Edges = {total_edges_from_corners} + {total_edges_from_edges} + {total_edges_from_interior} = {total_edges}")