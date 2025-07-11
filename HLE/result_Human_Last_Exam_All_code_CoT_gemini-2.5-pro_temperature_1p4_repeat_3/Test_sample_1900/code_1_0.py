# This script determines and prints the computational complexity for the two questions asked.

# The reasoning for the complexity is as follows:
# Let G = (V, E) be the graph where V is the set of n red balls and E
# represents the neighbor relationship.

# The properties of the red balls are translated into graph theory:
# 1. The graph G is connected.
# 2. For any vertex A, the subgraph induced by its neighbors, G[N(A)], is connected.
#    This defines a "neighbourhood-connected" graph.

# -- Analysis for Question A --
# Question A asks for the complexity of DECIDING if a path through all red balls exists.
# This is equivalent to deciding if the graph G has a Hamiltonian path.
# A mathematical theorem proves that any connected, neighbourhood-connected graph
# on n >= 2 vertices is guaranteed to have a Hamiltonian path.
# Since existence is guaranteed by the problem's premises, the decision is trivial
# and does not require any computation dependent on n.
# The complexity is constant time, O(1).

# -- Analysis for Question B --
# Question B asks for the complexity of FINDING such a path, given it exists.
# This requires an algorithm to construct a Hamiltonian path.
# The graph G is neighbourhood-connected, which is a subclass of "claw-free" graphs.
# Algorithms exist to find Hamiltonian paths in these special graph classes
# in polynomial time. A prominent complexity for such algorithms, for example,
# for claw-free graphs, is O(n^3).
# Therefore, the complexity for finding the path is O(n^3).

# Combine the answers into the specified format.
complexity_A = "O(1)"
complexity_B = "O(n^3)"

# Print the final answer.
print(f"{complexity_A}; {complexity_B}")