# This script calculates and prints the coefficients for the number of 
# closed tree-like walks of length 6 in a simple graph X.

# The expression is of the form:
# c_1 * e + c_2 * k + c_3 * p + c_4 * sum(deg(v) choose 2) + c_5 * sum(deg(v) choose 3)
# where e is the number of edges, k is the number of triangles, p is the number of P_4 paths.

# The coefficients are determined by counting the number of walks whose underlying 
# simple graph of edges is a tree. The total length of the walk is 6.

# c_1: Coefficient for e (P_2 subgraphs).
# The underlying tree has 1 edge, traversed 6 times.
# For each edge, there are 2 possible walks (one starting at each vertex).
c_1 = 2

# c_2: Coefficient for k (K_3 subgraphs).
# A K_3 (triangle) is not a tree, so it cannot be an underlying graph.
c_2 = 0

# c_3: Coefficient for p (P_4 subgraphs).
# The underlying tree is a P_4 (3 edges). Each edge is traversed twice.
# A combinatorial count reveals 6 such walks for each P_4 subgraph.
c_3 = 6

# c_4: Coefficient for sum(deg(v) choose 2) (P_3 subgraphs).
# The underlying tree is a P_3 (2 edges). One edge is traversed 4 times, the other twice.
# A combinatorial count reveals 12 such walks for each P_3 subgraph.
c_4 = 12

# c_5: Coefficient for sum(deg(v) choose 3) (K_{1,3} subgraphs).
# The underlying tree is a K_{1,3} (3 edges). Each edge is traversed twice.
# A combinatorial count reveals 12 such walks for each K_{1,3} subgraph.
c_5 = 12

# Printing the coefficients in the specified order.
print(f"c_1 = {c_1}")
print(f"c_2 = {c_2}")
print(f"c_3 = {c_3}")
print(f"c_4 = {c_4}")
print(f"c_5 = {c_5}")