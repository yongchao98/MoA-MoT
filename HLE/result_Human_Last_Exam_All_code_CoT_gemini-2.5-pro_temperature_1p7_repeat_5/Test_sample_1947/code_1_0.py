# Based on the combinatorial analysis, the coefficients are determined as follows:
# c1: From walks on a single edge (P_2). Walk is back-and-forth 3 times.
c1 = 2

# c2: From walks on a K_3. A K_3 is not a tree, so by definition, there are no tree-like walks whose edge-set forms a K_3.
c2 = 0

# c3: From walks on a P_4 subgraph. One round-trip on each of the 3 edges.
c3 = 6

# c4: From walks on a P_3 subgraph. One edge has 2 round-trips, the other has 1.
c4 = 10

# c5: From walks on a K_{1,3} subgraph. One round-trip on each of the 3 edges.
c5 = 12

# The problem asks for the coefficients c_1, c_2, c_3, c_4, c_5 in this order.
print(f"{c1},{c2},{c3},{c4},{c5}")