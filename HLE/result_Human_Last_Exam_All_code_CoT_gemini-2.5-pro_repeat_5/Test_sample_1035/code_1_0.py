import math

# The number of vertices in the complete graph K_n
n = 35

# Step 1: Find the minimum number of bipartite graphs (k) required to cover K_n.
# This is given by the formula k = ceil(log2(n)).
k = math.ceil(math.log2(n))

# Step 2: Calculate the total number of vertices in the family of graphs.
# The standard construction uses k spanning bipartite graphs, so each graph has n vertices.
# The total number of vertices is the sum of vertices in each graph, which is k * n.
total_vertices = k * n

# Print the final equation with each number.
print(f"The minimum number of bipartite graphs required is k = ceil(log2({n})) = {k}.")
print("The standard construction uses k spanning graphs.")
print(f"The minimum total number of vertices is k * n, which is:")
print(f"{k} * {n} = {total_vertices}")
