# Step 1: Define the folding parameters
num_top_bottom_folds = 2
num_left_right_folds = 2
total_folds = num_top_bottom_folds + num_left_right_folds

# Step 2: Calculate the edges from internal holes
# A cut at the center of the paper is unfolded by all 4 folds.
edges_center_hole = 2**total_folds

# A cut at the midpoints of the left/right edges lies on the horizontal centerline.
# It's unfolded by the 2 top-bottom folds, creating 2 holes.
edges_horiz_holes = 2 * (2**num_top_bottom_folds)

# A cut at the midpoints of the top/bottom edges lies on the vertical centerline.
# It's unfolded by the 2 left-right folds, creating 2 holes.
edges_vert_holes = 2 * (2**num_left_right_folds)

total_hole_edges = edges_center_hole + edges_horiz_holes + edges_vert_holes

# Step 3: Calculate the edges from the outer boundary
# The 4 original corners are cut, replacing each point with 1 edge.
edges_from_corners = 4

# The 2 midpoints on horizontal edges get notched. Each notch has 2^num_top_bottom_folds edges.
edges_from_horiz_notches = 2 * (2**num_top_bottom_folds)

# The 2 midpoints on vertical edges get notched. Each notch has 2^num_left_right_folds edges.
edges_from_vert_notches = 2 * (2**num_left_right_folds)

total_notch_edges = edges_from_horiz_notches + edges_from_vert_notches

# Each of the 4 original edges is broken into 2 segments between the corner cuts and midpoint notches.
edges_from_segments = 4 * 2

total_outer_edges = edges_from_corners + total_notch_edges + edges_from_segments

# Step 4: Calculate the total number of edges
total_edges = total_hole_edges + total_outer_edges

# Step 5: Print the final equation and result
print("Total Edges = (Central Hole Edges) + (Horizontal Holes Edges) + (Vertical Holes Edges) + (Outer Corner Edges) + (Outer Notch Edges) + (Outer Segment Edges)")
print(f"Total Edges = {edges_center_hole} + {edges_horiz_holes} + {edges_vert_holes} + {edges_from_corners} + ({edges_from_horiz_notches} + {edges_from_vert_notches}) + {edges_from_segments}")
print(f"Total Edges = {edges_center_hole} + {edges_horiz_holes} + {edges_vert_holes} + {edges_from_corners} + {total_notch_edges} + {edges_from_segments}")
print(f"Total Edges = {total_edges}")

print("<<<60>>>")