# The reasoning above leads to a specific number.
# The two possible graphs are determined by whether their "base graph" H on 1000 vertices
# is a cycle or a path.
# Case 1: H is a cycle (C_1000). This leads to one graph (the prism graph Y_1000).
# Case 2: H is a path (P_1000). This leads to one graph (a form of Möbius ladder graph).
# These two graphs are non-isomorphic.
# Total number of non-isomorphic graphs is 1 + 1 = 2.

num_graphs = 2
print(f"There are {num_graphs} non-isomorphic graphs with the given properties.")
print("The first type of graph is the prism graph Y_1000.")
print("The second type is constructed from two 1000-vertex paths, with their endpoints cross-connected.")
print(f"Final answer: {num_graphs}")