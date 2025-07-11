import math

# The orders of the fundamental groups of X1, X2, and X3
n1 = 5  # |pi_1(X1)| = |Z_5|
n2 = 8  # |pi_1(X2)| = |Z_8|
n3 = 2  # |pi_1(X3)| = |Z_2|

print(f"The orders of the factor groups are n1 = {n1}, n2 = {n2}, n3 = {n3}.")

# The number of vertices v in the quotient graph K \ T is the sum of the sizes
# of the quotients of the abelianization A by each factor group.
# v = |A/G1| + |A/G2| + |A/G3|
v = (n2 * n3) + (n1 * n3) + (n1 * n2)

print(f"The number of vertices in the quotient graph is v = ({n2}*{n3}) + ({n1}*{n3}) + ({n1}*{n2}) = {v}")

# The number of edges e in the quotient graph K \ T depends on the chosen
# graph of groups. For a line graph with 3 vertices, there are 2 edges.
# e = 2 * |A|
e = 2 * (n1 * n2 * n3)

print(f"The number of edges in the quotient graph is e = 2 * ({n1}*{n2}*{n3}) = {e}")

# The rank of the kernel K, as a free group, is the rank of the fundamental
# group of the quotient graph, which is rank = e - v + 1 (assuming it's connected).
rank = e - v + 1

print(f"The rank of the kernel K is given by the formula: Rank = e - v + 1.")
print(f"Rank = {e} - {v} + 1 = {rank}")
