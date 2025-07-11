#
# This script determines the computational complexity for the two questions posed.
#
# Problem Analysis:
# The problem of finding a "line through all the red balls" is equivalent to finding a
# Hamiltonian Path in a graph. The vertices of the graph are the red balls, and an
# edge exists between two vertices if the corresponding balls are neighbors.
# Let n be the number of red balls (vertices).
#
# The given properties of the set of red balls translate to:
# 1. The graph is connected.
# 2. The graph is locally connected (the subgraph induced by the neighbors of any vertex is connected).
#
# Question A: Complexity of deciding if a Hamiltonian Path exists.
# The Hamiltonian Path Problem is NP-complete, even for the special class of connected,
# locally connected graphs. This means there is no known algorithm that can solve this
# problem in polynomial time for all cases. The best-known exact algorithms, such as
# the dynamic programming approach by Held and Karp, have a worst-case time complexity
# that is exponential in the number of vertices, n.
#
# Question B: Complexity of finding a Hamiltonian Path.
# This is the search version of the problem. Finding the path is at least as hard as
# deciding if one exists. The same algorithm used to solve the decision problem can be
# adapted to reconstruct the path, and it does so within the same time complexity bound.
#
# The complexity for both problems is therefore dominated by an exponential term. The
# standard and most precise expression for the Held-Karp algorithm's complexity is
# O(n^2 * 2^n).

# The complexity for deciding if a path exists.
complexity_A = "O(n^2 * 2^n)"

# The complexity for finding the path.
complexity_B = "O(n^2 * 2^n)"

# Print the final answer in the format "O(...); O(...)".
# This output includes the numbers '2' as per the instructions.
print(f"{complexity_A}; {complexity_B}")
