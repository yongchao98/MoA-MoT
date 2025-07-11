# 1. Analyze the shapes created by the cuts and folds.
# The original paper is a square. After two pairs of folds (top-to-bottom, left-to-right),
# we have a small, 16-layer square. The four corners of this folded square are not
# equivalent in terms of their origin on the initial paper.
# - One corner corresponds to the four corners of the original square.
# - One corner corresponds to the exact center of the original square.
# - The other two corners correspond to points along the central fold lines (axes) of the original square.
#
# Cutting these four distinct corners results in:
# - A modification of the outer boundary of the paper.
# - One central hole.
# - Two sets of identical holes located symmetrically on the paper.

# 2. Calculate the number of edges on the outer boundary.
# The cut at the corner corresponding to the original square's corners effectively
# cuts off all four corners of the original paper.
# A square has 4 edges. Cutting off the 4 corners turns it into an octagon.
outer_boundary_edges = 8
print(f"The main shape's outer boundary results from cutting the four corners of the original square.")
print(f"This changes the 4-sided square into an 8-sided octagon.")
print(f"Number of edges on the outer boundary = {outer_boundary_edges}\n")

# 3. Calculate the number of edges for each hole.
# The number of edges of a hole created by cutting a folded corner is 2^k,
# where k is the number of folds that created that corner.

# 3a. The central hole.
# This is created by cutting the corner that corresponds to the paper's center.
# This point is affected by all 4 folds. So, k=4.
num_central_holes = 1
folds_for_central_hole = 4
central_hole_edges = 2**folds_for_central_hole
print(f"One hole is created at the very center of the paper.")
print(f"This corner was created by {folds_for_central_hole} folds, so the hole has 2^{folds_for_central_hole} edges.")
print(f"Edges from the central hole = {central_hole_edges}\n")

# 3b. The axis holes.
# The two remaining corners of the folded square lie on the main fold lines (axes).
# Each of these corners is created by 2 orthogonal folds (e.g., the first T-B fold and the fourth L-R fold).
# A cut at such a corner creates a hole with 2^k = 2^2 = 4 edges.
# Due to the symmetry of the other folds, the cut at each of these two distinct "axis" corners
# creates a set of 2 identical holes each.

# First set of axis holes
num_axis_holes_1 = 2
folds_for_axis_hole_1 = 2
edges_per_axis_hole_1 = 2**folds_for_axis_hole_1
total_axis_edges_1 = num_axis_holes_1 * edges_per_axis_hole_1
print(f"The cut on the third corner of the folded square creates a set of {num_axis_holes_1} holes.")
print(f"Each of these holes has 2^{folds_for_axis_hole_1} = {edges_per_axis_hole_1} edges.")
print(f"Total edges from this set of holes = {num_axis_holes_1} * {edges_per_axis_hole_1} = {total_axis_edges_1}\n")

# Second set of axis holes
num_axis_holes_2 = 2
folds_for_axis_hole_2 = 2
edges_per_axis_hole_2 = 2**folds_for_axis_hole_2
total_axis_edges_2 = num_axis_holes_2 * edges_per_axis_hole_2
print(f"The cut on the fourth corner creates another set of {num_axis_holes_2} holes.")
print(f"Each of these holes also has 2^{folds_for_axis_hole_2} = {edges_per_axis_hole_2} edges.")
print(f"Total edges from this set of holes = {num_axis_holes_2} * {edges_per_axis_hole_2} = {total_axis_edges_2}\n")


# 4. Sum all the edges.
total_edges = outer_boundary_edges + central_hole_edges + total_axis_edges_1 + total_axis_edges_2
print("The total number of edges is the sum of the edges from the outer boundary and all the holes.")
print(f"Final Equation: {outer_boundary_edges} + {central_hole_edges} + {total_axis_edges_1} + {total_axis_edges_2} = {total_edges}")