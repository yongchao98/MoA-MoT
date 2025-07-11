import math

# The problem states d is an even integer.
# The edge connectivity of G is 2, which implies the minimum degree of G is at least 2.
# So, d >= 2.
# We will use the smallest possible value for d, which is 2.
d = 2

# The total number of edges removed is the sum of the degrees of v1, v2, v3.
# deg(v1) = d
# deg(v2) = d + 1
# deg(v3) = d + 1
total_degree_sum = d + (d + 1) + (d + 1)
print(f"For d = {d}, the sum of degrees of the removed vertices is:")
print(f"{d} + {d + 1} + {d + 1} = {total_degree_sum}")

# The worst-case for G' is that it's maximally fragmented.
# Each fragment must have been connected to {v1, v2, v3} by at least 2 edges
# for the original graph G to be 2-edge-connected.
# So, the maximum number of leaf blocks (l) is total_degree_sum / 2.
l = total_degree_sum / 2
print(f"The maximum number of leaf blocks 'l' in G' is:")
print(f"{total_degree_sum} / 2 = {l}")

# The minimum number of edges to add to make a graph 2-edge-connected
# is ceil(l / 2).
num_edges_to_add = math.ceil(l / 2)
print(f"The minimal number of edges to add to G' is ceil(l / 2):")
print(f"ceil({l} / 2) = {num_edges_to_add}")
