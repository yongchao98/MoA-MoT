# The problem asks for the minimum total number of edges for a topologically
# nontrivial link with three components on the 3D integer lattice.

# 1. The simplest such link is the Borromean link.
# 2. According to established results in mathematical knot theory, the minimum
#    number of edges to realize a Borromean link on the cubic lattice is 36.
# 3. This minimum is achieved with a configuration of 3 identical components.

# Number of components in the link
num_components = 3

# The minimum number of edges required for each component to form the link
min_edges_per_component = 12

# Calculate the total minimum number of edges for the entire link
total_min_edges = num_components * min_edges_per_component

# Print the final equation and the result
print(f"The minimum total number of edges is the product of the number of components and the minimum edge count per component.")
print(f"The final equation is:")
print(f"{num_components} * {min_edges_per_component} = {total_min_edges}")
