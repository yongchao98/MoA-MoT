# The problem is to find the coefficients c_1, c_2, c_3, c_4, c_5 in the given expression
# for the number of closed tree-like walks of length 6.
# The coefficients are derived from a combinatorial analysis of the possible walk structures.

# c_1 is the coefficient for 'e', the number of edges.
# A walk over a single edge traversed 6 times. For each edge, there are 2 such walks.
c_1 = 2

# c_2 is the coefficient for 'k', the number of triangles (K_3).
# A walk whose edges form a K_3 is not tree-like, so they are not counted.
c_2 = 0

# c_3 is the coefficient for 'p', the number of P_4 subgraphs.
# A walk on a P_4 tree, with each of its 3 edges traversed twice. There are 6 such walks.
c_3 = 6

# c_4 is the coefficient for the sum of deg(v) choose 2, which counts P_3 subgraphs.
# A walk on a P_3 tree, with one edge traversed 4 times and the other twice. There are 12 such walks.
c_4 = 12

# c_5 is the coefficient for the sum of deg(v) choose 3, which counts K_{1,3} subgraphs.
# A walk on a K_{1,3} tree, with each of its 3 edges traversed twice. There are 12 such walks.
c_5 = 12

# Printing the coefficients in the specified order.
print(c_1)
print(c_2)
print(c_3)
print(c_4)
print(c_5)