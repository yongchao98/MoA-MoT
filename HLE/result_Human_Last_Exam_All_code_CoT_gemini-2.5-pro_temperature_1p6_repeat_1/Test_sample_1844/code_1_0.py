# Step 1: Define the number of edges for the outer boundary.
# A square has 4 edges. When one corner is cut off, two edges are shortened and one new edge is added,
# resulting in a pentagon shape for the outer boundary.
edges_outer_boundary = 5

# Step 2: Calculate the edges from the hole created by cutting the corner on one crease.
# This cut goes through 2^1 = 2 layers. Unfolding creates a single diamond-shaped hole.
# A diamond has 4 edges.
edges_hole_from_2_layers_1 = 4

# Step 3: Calculate the edges from the hole created by cutting the other corner on one crease.
# This cut also goes through 2^1 = 2 layers and creates another 4-edged hole.
edges_hole_from_2_layers_2 = 4

# Step 4: Calculate the edges from the hole created by cutting the corner on two creases.
# This cut goes through 2^2 = 4 layers. Unfolding across two perpendicular folds also creates
# a single diamond-shaped hole with 4 edges.
edges_hole_from_4_layers = 4

# Step 5: Sum all the edges to find the total.
total_edges = edges_outer_boundary + edges_hole_from_2_layers_1 + edges_hole_from_2_layers_2 + edges_hole_from_4_layers

# Output the breakdown of the calculation.
print("The final shape has an outer boundary and three internal holes.")
print(f"The outer boundary has {edges_outer_boundary} edges.")
print(f"The first hole has {edges_hole_from_2_layers_1} edges.")
print(f"The second hole has {edges_hole_from_2_layers_2} edges.")
print(f"The third hole has {edges_hole_from_4_layers} edges.")
print(f"The total number of edges is the sum: {edges_outer_boundary} + {edges_hole_from_2_layers_1} + {edges_hole_from_2_layers_2} + {edges_hole_from_4_layers} = {total_edges}")