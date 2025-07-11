# The number of vertices in the graph G
n = 12345

# A graph on n vertices has a chromatic number of n if and only if it is the complete graph Kn.
# The problem states that G is not the complete graph.
# Therefore, the chromatic number of G must be less than n.
# The maximum integer value less than n is n - 1. So, an upper bound for the chromatic number is n - 1.

# To verify this is the maximum possible value, we need to show that there exists a non-complete graph
# on n vertices that requires n - 1 colors.
# Consider the graph formed by taking the complete graph Kn and removing a single edge.
# This graph is not complete, and it can be shown that its chromatic number is exactly n - 1.
# (As explained in the steps above, you can color the two non-adjacent vertices with the same color,
# and all other n-2 vertices will need their own unique colors).

# So, the maximum number of colours is n - 1.
max_colors = n - 1

# Print the equation and the final answer.
print(f"The maximum number of colors is found by the equation:")
print(f"{n} - 1 = {max_colors}")