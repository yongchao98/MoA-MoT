def solve_clique_packing():
    """
    Calculates the maximum number of different-sized, vertex-disjoint cliques
    that can fit into a graph with a given number of vertices.
    """
    n = 128  # Total number of vertices in the graph

    # We want to find the maximum number of different clique sizes, let's call it 'm'.
    # To maximize 'm', we should choose the smallest possible distinct clique sizes: 1, 2, 3, ...
    # We need to find the largest 'm' such that 1 + 2 + ... + m <= n.

    clique_sizes = []
    total_vertices_used = 0
    current_clique_size = 1

    # Greedily add cliques of increasing size as long as they fit.
    while total_vertices_used + current_clique_size <= n:
        # "Pack" a clique of the current size
        total_vertices_used += current_clique_size
        clique_sizes.append(current_clique_size)

        # Move to the next distinct size
        current_clique_size += 1

    # The maximum number of different clique sizes is the number of terms we could sum up.
    max_m = len(clique_sizes)

    # As requested, we will output the final equation.
    # This equation shows the sum of the sizes of the cliques we packed.
    equation_str = " + ".join(map(str, clique_sizes))
    final_sum = sum(clique_sizes)

    print(f"To find the maximum number of different clique sizes (m) that can appear as vertex-disjoint induced subgraphs in a graph with {n} vertices, we find the largest m such that the sum of the first m positive integers is at most {n}.")
    print("\nThe final equation for the sum of the chosen clique sizes is:")
    print(f"{equation_str} = {final_sum}")
    print(f"\nThe total vertices used is {final_sum}, which is <= {n}.")
    print(f"The next clique size would be {current_clique_size}, making the sum {final_sum + current_clique_size}, which exceeds {n}.")
    print(f"\nThus, the maximum possible number of different clique sizes is {max_m}.")


solve_clique_packing()

# The final answer is the maximum number of different clique sizes.
# Based on the calculation, this is 15.
# 1 + 2 + ... + 15 = 120 <= 128
# 1 + 2 + ... + 16 = 136 > 128
print("\n<<<15>>>")