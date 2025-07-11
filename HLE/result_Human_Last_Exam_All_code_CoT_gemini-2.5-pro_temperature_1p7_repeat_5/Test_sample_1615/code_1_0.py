# The number of vertices in the graph.
n = 12345

# For any graph G with n vertices, its chromatic number chi(G) is at most n.
# The only graph requiring n colors is the complete graph K_n.
# The problem states that G is not a complete graph, so its chromatic number must be less than n.
# This means the maximum possible chromatic number is at most n - 1.

# To show this maximum can be reached, consider a graph G' created by removing a single edge from K_n.
# This graph G' is not complete, and its chromatic number is exactly n - 1.
# Therefore, the maximum number of colors that might be needed is n - 1.

max_colors = n - 1

# We are asked to output the final equation.
print("The maximum number of colours is given by the equation:")
print(f"{n} - 1 = {max_colors}")