# This script calculates the total number of edges on the paper after unfolding.

# Part 1: Calculate edges from the internal hole.
# The cut at the corner formed by two fold lines is at the center of the original paper.
# This creates a single internal hole. The paper is folded 4 times. An internal cut's
# edge count is doubled with each unfolding. The initial cut creates 1 edge.
num_folds = 4
initial_cut_edge = 1
edges_hole = initial_cut_edge * (2**num_folds)

# Part 2: Calculate edges on the outer boundary.
# We start with an uncut square, which has 4 edges.
base_outer_edges = 4

# The cut at the original paper's corner is not on a fold line, so it's a single
# snip that replaces a point (0 edges) with a line (1 edge).
edge_add_from_corner_cut = 1

# The other two cuts are on the edges of the folded square. Each corresponds to a
# point on an edge of the original paper (e.g., one on the bottom edge, one on the right edge).
# - A cut on a folded edge is mirrored upon unfolding, creating a V-shaped indentation with 2 edges.
# - This indentation is then replicated onto the opposite side of the paper by another fold.
# - This results in 2 V-shaped indentations from each of these two initial cuts (4 total indentations).
# A V-indentation replaces 1 straight edge segment with 4 edges (2 for the 'V' and 2 remaining parts
# of the original edge). This is a net gain of 3 edges per indentation.
num_indentations = 4
net_edge_add_per_indentation = 3
edge_add_from_indentations = num_indentations * net_edge_add_per_indentation

# The total number of edges on the outer boundary is the sum of the initial edges and the net additions.
edges_outer = base_outer_edges + edge_add_from_corner_cut + edge_add_from_indentations

# Part 3: Calculate the total number of edges.
total_edges = edges_hole + edges_outer

# Print the breakdown of the calculation as requested.
print("The total number of edges is the sum of edges from the internal hole and the outer boundary.")
print(f"Number of edges from the central hole: {edges_hole}")
print(f"Number of edges on the outer boundary: {edges_outer}")
print("The final calculation is:")
print(f"{edges_hole} + {edges_outer} = {total_edges}")
