import math

# The total number of vertices available in the graph.
n = 128

# Variables to keep track of the calculation.
total_vertices_used = 0
num_different_sizes = 0
clique_sizes = []

# This loop determines the maximum number of different clique sizes (m)
# by summing the smallest possible distinct sizes (1, 2, 3, ...)
# until the total number of vertices required exceeds n.
while True:
    # The next smallest size to try.
    next_size = num_different_sizes + 1
    
    # Check if adding a clique of this new size is possible.
    if total_vertices_used + next_size <= n:
        total_vertices_used += next_size
        clique_sizes.append(next_size)
        num_different_sizes += 1
    else:
        # If we cannot add a clique of this size, we have found our maximum.
        break

print(f"The task is to find the maximum number of different maximal clique sizes in a graph with n = {n} vertices.")
print("To achieve this, we construct a graph composed of disjoint cliques of the smallest possible different sizes (1, 2, 3, ...).")
print("We need to find the largest number of such cliques whose total vertices do not exceed 128.")
print("\nThe calculation involves summing the clique sizes until the limit is reached.")

# Format the equation string as requested.
equation_str = " + ".join(map(str, clique_sizes))

print("\nThe final equation representing the sum of vertices for the chosen cliques is:")
print(f"{equation_str} = {total_vertices_used}")

remaining_vertices = n - total_vertices_used
print(f"\nThis construction uses {total_vertices_used} vertices out of {n}, leaving {remaining_vertices} unused.")
print(f"The next clique size to consider would be {num_different_sizes + 1}.")
print(f"Adding it would require a total of {total_vertices_used} + {num_different_sizes + 1} = {total_vertices_used + num_different_sizes + 1} vertices, which exceeds {n}.")
print(f"\nTherefore, the maximum possible number of different clique sizes is the number of terms in the sum.")
print(f"The answer is {num_different_sizes}.")
