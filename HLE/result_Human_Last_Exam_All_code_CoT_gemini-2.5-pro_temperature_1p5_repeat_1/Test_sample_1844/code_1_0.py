# 1. Calculate the number of internal edges from the "center" cut.
# This cut creates internal holes when unfolded.
num_center_corners = 1
# This single cut results in 4 separate diamond-shaped holes.
holes_per_center_cut = 4
# Each diamond-shaped hole has 4 edges.
edges_per_hole = 4
internal_edges = num_center_corners * holes_per_center_cut * edges_per_hole

print(f"Calculation for internal edges (holes):")
print(f"The cut near the 'center' corner creates {holes_per_center_cut} separate holes.")
print(f"Each hole is a diamond shape and has {edges_per_hole} edges.")
print(f"Total internal edges = {num_center_corners} cut * {holes_per_center_cut} holes/cut * {edges_per_hole} edges/hole = {internal_edges}")
print("-" * 20)

# 2. Calculate the number of edges on the outer boundary.
# The original square has 4 edges.
initial_outer_edges = 4
print(f"Calculation for the outer boundary edges:")
print(f"The original square has {initial_outer_edges} edges.")

# The "corner" cut unfolds to clip the 4 corners of the square, making an octagon.
# An octagon has 8 edges. This adds 4 new edges to the boundary count.
edges_added_by_corner_clips = 4
edges_after_clipping = initial_outer_edges + edges_added_by_corner_clips
print(f"Clipping the 4 corners adds {edges_added_by_corner_clips} edges, turning the square into an octagon with {edges_after_clipping} edges.")

# The two "edge" cuts unfold into 4 notches on the sides of the shape.
# Each notch adds 1 edge to the total boundary count.
num_edge_corners = 2
notches_per_edge_corner_cut = 2
edge_increase_per_notch = 1
total_edges_added_by_notches = num_edge_corners * notches_per_edge_corner_cut * edge_increase_per_notch
outer_edges = edges_after_clipping + total_edges_added_by_notches
print(f"The 2 'edge' cuts create {num_edge_corners * notches_per_edge_corner_cut} notches on the boundary.")
print(f"Each notch adds {edge_increase_per_notch} edge to the count.")
print(f"Total outer edges = {edges_after_clipping} (octagon) + {total_edges_added_by_notches} (from notches) = {outer_edges}")
print("-" * 20)

# 3. Sum internal and outer edges for the total.
total_edges = internal_edges + outer_edges
print("Final Calculation:")
print(f"Total Edges = Internal Edges + Outer Edges")
print(f"Total Edges = {internal_edges} + {outer_edges} = {total_edges}")
print("-" * 20)
print(f"The total number of edges on the final shape is {total_edges}.")
