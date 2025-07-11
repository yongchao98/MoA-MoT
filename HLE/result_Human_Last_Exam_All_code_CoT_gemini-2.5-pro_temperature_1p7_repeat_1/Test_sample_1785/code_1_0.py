# The problem asks for the minimum total number of edges in a topologically
# nontrivial link with three components on a 3D integer lattice.
# This corresponds to the minimal lattice realization of the Borromean rings.

# From established results in knot theory, it is known that the minimal
# representation of the Borromean rings on a cubic lattice consists of
# three congruent, non-planar components.

# Each of these minimal components requires 12 edges.
edges_K1 = 12
edges_K2 = 12
edges_K3 = 12

# The total number of edges is the sum of the edges in the three components.
total_edges = edges_K1 + edges_K2 + edges_K3

# The final output prints the calculation, showing the number of edges for each
# component and the final sum, as requested.
print(f"The minimum total number of edges is found by summing the edges of the three minimal components:")
print(f"{edges_K1} + {edges_K2} + {edges_K3} = {total_edges}")
