def solve_treewidth_bound():
    """
    Calculates the tight upper bound on the treewidth of the union of two graphs.
    """
    # Example values for the treewidths of H and G, and the number of shared vertices.
    # The bound is independent of k, but we include it to reflect the problem statement.
    t_H = 12
    t_G = 9
    k = 5

    # The tight upper bound is max(t_H, t_G).
    bound = max(t_H, t_G)

    # Output the explanation and the result.
    print(f"Given the treewidth of graph H, t_H = {t_H}")
    print(f"Given the treewidth of graph G, t_G = {t_G}")
    print(f"Given the number of shared vertices, k = {k}")
    print("\nThe tight upper bound on the treewidth of the combined graph F = H U G is given by max(t_H, t_G).")
    print("\nPlugging in the values, we get the final equation:")
    # Print the final equation with numbers as requested
    print(f"max({t_H}, {t_G}) = {bound}")

solve_treewidth_bound()