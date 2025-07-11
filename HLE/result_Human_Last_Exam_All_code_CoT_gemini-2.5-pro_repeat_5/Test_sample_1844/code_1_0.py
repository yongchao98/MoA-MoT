# This script calculates the total number of edges on the paper after it's been folded, cut, and unfolded.

# Step 1: Analyze the cuts based on their location on the original paper.
# The paper is folded four times, creating 2^4 = 16 layers.
# The four corners of the folded square correspond to:
# 1. The absolute center of the original paper.
# 2. The four original corners of the paper, stacked together.
# 3. The midpoints of the top and bottom original edges, stacked together.
# 4. The midpoints of the left and right original edges, stacked together.

print("Calculating the total number of edges for the unfolded shape.")
print("----------------------------------------------------------")

# Step 2: Calculate the edges of the internal hole.
# The cut at the corner corresponding to the paper's center creates a hole.
# When unfolded, a single straight cut at the center is mirrored across two perpendicular fold lines,
# forming a closed shape with four sides.
hole_edges = 4
print(f"The cut at the center creates 1 internal hole.")
print(f"The number of edges for the hole is: {hole_edges}")
print(" ")

# Step 3: Calculate the edges on the outer perimeter.
print("Calculating the edges on the outer perimeter...")

# a) Edges from cutting the 4 original corners.
# The cut at the corner corresponding to the original four corners clips all of them simultaneously.
# Each clipped corner replaces a point with one new edge.
num_original_corners = 4
edges_from_corner_cuts = 1 * num_original_corners
print(f" - Each of the {num_original_corners} original corners is cut, creating 1 new edge per corner.")
print(f"   Edges from corner cuts = {edges_from_corner_cuts}")

# b) Edges from cutting notches in the middle of the 4 original sides.
# The two cuts at the 'midpoint' corners create V-shaped notches at the midpoint of each of the 4 original sides.
# Each V-shaped notch consists of 2 edges.
num_side_notches = 4
edges_per_notch = 2
edges_from_notches = num_side_notches * edges_per_notch
print(f" - Each of the {num_side_notches} original sides gets a notch, creating {edges_per_notch} new edges per notch.")
print(f"   Edges from notches = {edges_from_notches}")

# c) Remaining segments of the original 4 edges.
# Each of the 4 original edges is split into 2 segments by the notch in its middle.
num_original_sides = 4
segments_per_side = 2
original_edge_segments = num_original_sides * segments_per_side
print(f" - Each of the {num_original_sides} original sides is split into {segments_per_side} segments.")
print(f"   Remaining original edge segments = {original_edge_segments}")
print(" ")

# Step 4: Calculate the total number of edges.
total_perimeter_edges = edges_from_corner_cuts + edges_from_notches + original_edge_segments
print(f"Total edges on the outer perimeter = {edges_from_corner_cuts} + {edges_from_notches} + {original_edge_segments} = {total_perimeter_edges}")

total_edges = total_perimeter_edges + hole_edges
print(f"Total edges (Perimeter + Hole) = {total_perimeter_edges} + {hole_edges}")
print("----------------------------------------------------------")
print(f"The final shape will have a total of {total_edges} edges.")
