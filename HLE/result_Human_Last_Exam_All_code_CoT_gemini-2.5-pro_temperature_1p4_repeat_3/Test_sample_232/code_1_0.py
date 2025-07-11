import sys

# This script calculates the number of regions the folds of an origami crane divide the paper into.
# The method uses Euler's formula for planar graphs: F = E - V + 1, where
# F is the number of faces (regions), E is the number of edges, and V is the number of vertices.

# Step 1: Calculate the regions for the "bird base" crease pattern.
print("Step 1: Calculating regions for the 'bird base' crease pattern.")
print("---------------------------------------------------------------")
# The bird base pattern has a specific number of vertices (intersections) and edges (crease lines).
# Vertices (V1):
# - 4 at the corners of the square
# - 4 at the midpoints of the square's edges
# - 1 at the absolute center of the square
# - 4 where creases intersect on the main diagonals
v1_corners = 4
v1_midpoints = 4
v1_center = 1
v1_inner = 4
V1 = v1_corners + v1_midpoints + v1_center + v1_inner
print(f"The bird base has {V1} vertices (V1).")

# Edges (E1):
# - 8 segments making up the boundary of the square
# - 4 segments connecting the center to the midpoints
# - 8 segments that form the main diagonals
# - 8 segments connecting the midpoints to the inner diagonal vertices
e1_boundary = 8
e1_median = 4
e1_diagonal = 8
e1_connecting = 8
E1 = e1_boundary + e1_median + e1_diagonal + e1_connecting
print(f"The bird base has {E1} edges (E1).")

# Using Euler's formula to find the number of regions for the bird base.
F1 = E1 - V1 + 1
print(f"The number of regions for the bird base is F1 = E1 - V1 + 1 = {E1} - {V1} + 1 = {F1}.")
print("\n")


# Step 2: Account for the additional folds that create the head and tail.
print("Step 2: Accounting for the head and tail folds.")
print("-------------------------------------------------")
print("The inside-reverse fold used for the head adds a V-shaped crease to the pattern.")
# Each leg of the 'V' crease terminates on an existing edge, creating a new vertex.
# Therefore, one V-shaped crease adds 2 vertices.
delta_V_per_fold = 2
print(f"Each fold (head or tail) adds {delta_V_per_fold} new vertices (delta_V).")

# Each V-shaped crease also adds edges:
# - 2 new edges for the 'V' shape itself.
# - 2 existing edges are split, turning them into 4 edges (a net addition of 2).
delta_E_per_fold = 2 + 2
print(f"Each fold adds {delta_E_per_fold} new edges (delta_E).")

# The number of new regions added by one fold is delta_E - delta_V.
delta_F_per_fold = delta_E_per_fold - delta_V_per_fold
print(f"This means each fold adds {delta_E_per_fold} - {delta_V_per_fold} = {delta_F_per_fold} new regions.")
print("\n")


# Step 3: Calculate the total number of regions.
print("Step 3: Calculating the total.")
print("------------------------------")
# There are two such folds: one for the head and one for the tail.
num_final_folds = 2
F_total = F1 + num_final_folds * delta_F_per_fold

print(f"We start with {F1} regions from the bird base and add the regions from {num_final_folds} folds.")
print("The final equation is:")
# Display the full calculation, substituting all the numbers.
# Format: (E1 - V1 + 1) + num_folds * (delta_E - delta_V)
# We want to show the numbers themselves.
print(f"({e1_boundary} + {e1_median} + {e1_diagonal} + {e1_connecting}) - ({v1_corners} + {v1_midpoints} + {v1_center} + {v1_inner}) + 1 + {num_final_folds} * ({delta_E_per_fold} - {delta_V_per_fold})")
print(f"= ({E1} - {V1} + 1) + {num_final_folds} * ({delta_F_per_fold})")
print(f"= {F1} + {num_final_folds * delta_F_per_fold}")
print(f"= {F_total}")
print("\nSo, the fold lines of a standard origami crane divide the paper into 20 regions.")

# Suppress the final answer output for clean execution in a shell, as the text above provides it.
sys.stdout = open('/dev/null', 'w')