# The coefficients for the expression for the number of closed tree-like walks of length 6.
# c1 for e (number of edges)
# c2 for k (number of K_3 subgraphs)
# c3 for p (number of P_4 subgraphs)
# c4 for sum over v of deg(v) choose 2 (number of P_3 subgraphs)
# c5 for sum over v of deg(v) choose 3 (number of K_1,3 subgraphs)

c1 = 2
c2 = 0
c3 = 6
c4 = 12
c5 = 6

# Print the coefficients in order.
print(f"{c1},{c2},{c3},{c4},{c5}")