# The problem is to find the minimum total number of edges for a topologically
# nontrivial link of three components on a 3D integer lattice. This specific
# type of link is known as the Borromean rings.

# Number of components in the link.
num_components = 3

# It is a known result from geometric knot theory that the minimum number of
# edges required for a single component to form a Borromean rings link on a
# lattice is 12. The minimal link is composed of three such identical components.
edges_per_component = 12

# Calculate the total minimum number of edges.
total_edges = num_components * edges_per_component

# Print the equation showing how the total is calculated, as requested.
print(f"The minimum total number of edges is found by multiplying the number of components by the minimum edges per component.")
print(f"{num_components} * {edges_per_component} = {total_edges}")