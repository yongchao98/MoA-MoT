import numpy as np

def get_trace_A6(adj_matrix):
    """Computes the trace of A^6 for a given adjacency matrix A."""
    A = np.array(adj_matrix, dtype=object)
    A6 = np.linalg.matrix_power(A, 6)
    return np.trace(A6)

# Case 1: T = P2 (a single edge)
# The graph has 2 vertices, 1 edge. Let's call it e1.
A_P2 = [[0, 1], [1, 0]]
C_P2 = get_trace_A6(A_P2)
# W(P2) = C(P2) since it has only one edge.
W_P2 = C_P2
c1 = W_P2
# Contribution is W(P2) * (number of P2 subgraphs) = c1 * e

# A K3 is not a tree. Its tree subgraphs are P3s and P2s.
# The contributions from these are already counted by other terms.
# So, the coefficient for k (number of K3) is 0.
c2 = 0

# Case 2: T = P3 (path with 2 edges)
# Vertices 0-1-2. Edges e1={0,1}, e2={1,2}.
A_P3 = [[0, 1, 0], [1, 0, 1], [0, 1, 0]]
C_P3 = get_trace_A6(A_P3)
# Subgraphs of P3 are P3 itself, two P2s (e1, e2), and the empty graph.
# W(P3) = C(P3) - (C(e1) + C(e2)) + C(empty)
# C(empty) is 0. C(e1)=C(e2)=C_P2.
W_P3 = C_P3 - 2 * C_P2
c4 = W_P3
# Contribution is W(P3) * (number of P3 subgraphs).
# Number of P3s is sum(deg(v) choose 2) for all v. So this is c4.

# Case 3: T = P4 (path with 3 edges)
# Vertices 0-1-2-3. Edges e1, e2, e3.
A_P4 = [[0, 1, 0, 0], [1, 0, 1, 0], [0, 1, 0, 1], [0, 0, 1, 0]]
C_P4 = get_trace_A6(A_P4)
# Subgraphs of P4 with edge sets: {e1,e2,e3}, {e1,e2}, {e2,e3}, {e1,e3}, {e1}, {e2}, {e3}, empty.
# C({e1,e2}) = C({e2,e3}) = C_P3.
# C({e1,e3}) is for two disjoint edges. C(G1 U G2) = C(G1)+C(G2)
C_P2_disjoint_P2 = C_P2 + C_P2
# W(P4) = C(P4) - (C({e1,e2}) + C({e2,e3}) + C({e1,e3})) + (C(e1)+C(e2)+C(e3)) - C(empty)
W_P4 = C_P4 - (C_P3 + C_P3 + C_P2_disjoint_P2) + (C_P2 + C_P2 + C_P2) - 0
c3 = W_P4
# Contribution is W(P4) * (number of P4 subgraphs) = c3 * p.

# Case 4: T = K1,3 (star graph with 3 edges)
# Center vertex 0, leaves 1,2,3. Edges e1,e2,e3.
A_K13 = [[0, 1, 1, 1], [1, 0, 0, 0], [1, 0, 0, 0], [1, 0, 0, 0]]
C_K13 = get_trace_A6(A_K13)
# Any pair of edges from K1,3 forms a P3. There are 3 such pairs.
# W(K1,3) = C(K1,3) - 3*C(P3) + 3*C(P2) - C(empty)
W_K13 = C_K13 - 3 * C_P3 + 3 * C_P2
c5 = W_K13
# Contribution is W(K1,3) * (number of K1,3 subgraphs).
# Number of K1,3s is sum(deg(v) choose 3) for all v. So this is c5.


# Print the final result in the required order.
# c1 for e, c2 for k, c3 for p, c4 for sum(deg choose 2), c5 for sum(deg choose 3).
print(f"c_1 = {c1}")
print(f"c_2 = {c2}")
print(f"c_3 = {c3}")
print(f"c_4 = {c4}")
print(f"c_5 = {c5}")

print("\nThe final coefficients c_1, c_2, c_3, c_4, c_5 are:")
print(f"{int(c1)}, {int(c2)}, {int(c3)}, {int(c4)}, {int(c5)}")