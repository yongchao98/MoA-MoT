# The problem asks for the minimum total number of edges in a topologically nontrivial
# 3-component link on a 3D integer lattice. This corresponds to finding the minimal
# lattice realization of the Borromean rings.

# According to established results in mathematical knot theory, the minimal
# configuration for the Borromean rings on a simple cubic lattice consists of
# three identical (congruent) components.

# The minimum number of edges required for each of these components is 12.
edges_per_component_1 = 12
edges_per_component_2 = 12
edges_per_component_3 = 12

# The total number of edges is the sum of the edges of the three components.
total_edges = edges_per_component_1 + edges_per_component_2 + edges_per_component_3

# We print the result as a clear equation.
print("The minimum total number of edges is the sum of the lengths of the three components.")
print(f"Component 1 Edges: {edges_per_component_1}")
print(f"Component 2 Edges: {edges_per_component_2}")
print(f"Component 3 Edges: {edges_per_component_3}")
print(f"Total Edges = {edges_per_component_1} + {edges_per_component_2} + {edges_per_component_3} = {total_edges}")
