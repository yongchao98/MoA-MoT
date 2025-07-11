# Based on the combinatorial analysis, the coefficients are determined as follows:
# c1: From walks on a single edge (P_2 tree).
# c2: From walks on K_3 subgraphs (no contribution as they are not tree-like).
# c3: From walks on P_4 subgraphs.
# c4: From walks on P_3 subgraphs.
# c5: From walks on K_{1,3} subgraphs.

c_1 = 2
c_2 = 0
c_3 = 12
c_4 = 12
c_5 = 12

# The expression for the number of closed tree-like walks of length 6 is:
# 2*e + 0*k + 12*p + 12*sum(deg(v) choose 2) + 12*sum(deg(v) choose 3)
# The problem asks for the coefficients in order.
print(c_1)
print(c_2)
print(c_3)
print(c_4)
print(c_5)