# Step 1: Define the initial state of the paper.
# A square has 4 edges.
original_outer_edges = 4

# Step 2: Calculate the edges of the outer boundary after the cuts.
# Of the four cuts on the folded paper, only one is at a "true" corner of the
# original paper. Cutting a corner of a polygon increases its edge count by one.
# For example, cutting a corner of a square (4 edges) makes it a pentagon (5 edges).
boundary_cut_edges = 1
final_outer_edges = original_outer_edges + boundary_cut_edges

# Step 3: Calculate the edges from the internal holes created by the other three cuts.
# The analysis of unfolding shows that the remaining three cuts create distinct sets of holes.

# The cut at the corner formed by an original edge and the 3rd fold line
# creates one 4-sided hole.
hole_set_1_count = 1
hole_set_1_sides = 4
hole_set_1_edges = hole_set_1_count * hole_set_1_sides

# The cut at the corner formed by an original edge and the 4th fold line
# creates two 4-sided holes.
hole_set_2_count = 2
hole_set_2_sides = 4
hole_set_2_edges = hole_set_2_count * hole_set_2_sides

# The cut at the corner formed by the intersection of the 3rd and 4th fold lines
# creates four 4-sided holes.
hole_set_3_count = 4
hole_set_3_sides = 4
hole_set_3_edges = hole_set_3_count * hole_set_3_sides

# Step 4: Calculate the total number of edges.
# This is the sum of the outer boundary edges and all internal hole edges.
total_internal_edges = hole_set_1_edges + hole_set_2_edges + hole_set_3_edges
total_edges = final_outer_edges + total_internal_edges

# Step 5: Print the final calculation step-by-step.
print(f"The final shape has an outer boundary with {final_outer_edges} edges.")
print(f"The first set of internal holes has {hole_set_1_edges} edges.")
print(f"The second set of internal holes has {hole_set_2_edges} edges.")
print(f"The third set of internal holes has {hole_set_3_edges} edges.")
print(f"The total number of edges is the sum of all these parts.")
print(f"Total Edges = {final_outer_edges} + {hole_set_1_edges} + {hole_set_2_edges} + {hole_set_3_edges} = {total_edges}")
