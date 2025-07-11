# Based on the step-by-step combinatorial analysis, the coefficients are determined as follows:
# c1: Contribution from walks on P_2 subgraphs (edges).
# c2: Contribution from walks related to K_3 subgraphs (triangles).
# c3: Contribution from walks on P_4 subgraphs (paths of length 3).
# c4: Contribution from walks on P_3 subgraphs (paths of length 2).
# c5: Contribution from walks on K_1,3 subgraphs (claws).

# Number of walks for each underlying tree structure:
# - On an edge (P_2), with 6 traversals: 2 walks.
# - On a P_3, with (4,2) and (2,4) traversals: 12 walks.
# - On a P_4, with (2,2,2) traversals: 6 walks.
# - On a K_1,3, with (2,2,2) traversals: 12 walks.

# The number of subgraphs corresponding to each coefficient are:
# - e (number of P_2)
# - k (number of K_3)
# - p (number of P_4)
# - sum(deg(v) choose 2) (number of P_3)
# - sum(deg(v) choose 3) (number of K_1,3)

# Thus, the coefficients are the number of walks per subgraph.
c1 = 2
c2 = 0  # A K_3 (triangle) is not a tree, so it cannot be the basis for a tree-like walk.
c3 = 6
c4 = 12
c5 = 12

# The problem requires printing the coefficients in order c_1, c_2, c_3, c_4, c_5.
print(f"{c1},{c2},{c3},{c4},{c5}")