# The number of vertices in the graph.
n = 128

# The problem asks for the maximum possible number of different clique sizes
# that can simultaneously appear as induced subgraphs in a graph with n vertices.

# The optimal construction is a complete graph, K_n.
# In a complete graph, any subset of k vertices forms an induced K_k clique.
# This means a K_n graph contains induced cliques of all possible integer sizes from 1 to n.
# The set of different clique sizes is {1, 2, 3, ..., n}.

# The total count of these distinct sizes is simply n.
maximum_number_of_sizes = n

# We can't have more than n sizes, because any clique's size k must be
# between 1 and n. The K_n construction shows that n is achievable.

# The final equation is that the maximum number of sizes is equal to n.
print(f"For a graph with n = {n} vertices:")
print(f"The maximum possible number of different clique sizes is {maximum_number_of_sizes}.")
print("The final equation is:")
print(f"{maximum_number_of_sizes} = {n}")
