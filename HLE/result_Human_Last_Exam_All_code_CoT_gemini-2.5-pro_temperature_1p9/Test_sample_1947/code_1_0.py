# The plan is to determine the coefficients by combinatorial counting
# for each type of tree-like walk, and then match them to the given expression.
# The thinking process is outlined above. The final coefficients derived are c1, c2, c3, c4, c5.

# c_1 corresponds to walks on a P_2 (single edge) subgraph.
# Contribution is 2 * e. So, c_1 = 2.
c1 = 2

# c_2 corresponds to the number of K_3 (triangles).
# A tree-like walk's underlying structure is a tree. The existence of an edge
# that would complete a triangle does not affect the count of walks on the
# tree structure itself. Therefore, the number of triangles does not directly
# contribute to the count of *tree-like* walks. So, c_2 = 0.
c2 = 0

# c_3 corresponds to walks on a P_4 (path of length 3) subgraph.
# Our combinatorial count yielded 6 walks per P_4.
# Contribution is 6 * p. So, c_3 = 6.
c3 = 6

# c_4 corresponds to walks on a P_3 (path of length 2) subgraph.
# Our combinatorial count yielded 12 walks per P_3.
# The number of P_3 subgraphs is given by sum(deg(v) choose 2).
# Contribution is 12 * sum(deg(v) choose 2). So, c_4 = 12.
c4 = 12

# c_5 corresponds to walks on a K_{1,3} (star) subgraph.
# Our combinatorial count yielded 12 walks per K_{1,3}.
# The number of K_{1,3} subgraphs is given by sum(deg(v) choose 3).
# Contribution is 12 * sum(deg(v) choose 3). So, c_5 = 12.
c5 = 12

print(f"{c1},{c2},{c3},{c4},{c5}")
