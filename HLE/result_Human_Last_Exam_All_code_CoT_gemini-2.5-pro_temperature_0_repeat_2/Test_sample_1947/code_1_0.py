# The coefficients for the expression for the number of closed tree-like walks of length 6.
# The expression is of the form:
# c_1 * e + c_2 * k + c_3 * p + c_4 * sum(deg(v) choose 2) + c_5 * sum(deg(v) choose 3)

# c_1: Coefficient for the number of edges (e)
# This corresponds to walks on a P_2 subgraph (a single edge) traversed 6 times.
# For each edge (u,v), there are two walks: u->v->u... and v->u->v...
c_1 = 2

# c_2: Coefficient for the number of K_3 subgraphs (k)
# A tree-like walk cannot have an underlying graph with a cycle. K_3 is a cycle.
# Thus, this term does not contribute.
c_2 = 0

# c_3: Coefficient for the number of P_4 subgraphs (p)
# This corresponds to walks on a P_4 subgraph, with each of the 3 edges traversed twice.
# Counting these walks gives 6 distinct walks for each P_4.
c_3 = 6

# c_4: Coefficient for the number of P_3 subgraphs (sum(deg(v) choose 2))
# This corresponds to walks on a P_3 subgraph. One edge is traversed 4 times, the other twice.
# For each P_3, there are two scenarios (which edge is traversed 4 times).
# Each scenario gives 6 walks, for a total of 12 walks per P_3.
c_4 = 12

# c_5: Coefficient for the number of K_{1,3} subgraphs (sum(deg(v) choose 3))
# This corresponds to walks on a K_{1,3} subgraph, with each of the 3 edges traversed twice.
# Counting these walks gives 12 distinct walks for each K_{1,3}.
c_5 = 12

# Print the coefficients in the specified order.
print(c_1, c_2, c_3, c_4, c_5)