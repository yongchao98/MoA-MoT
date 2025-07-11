def solve_treewidth_bound():
    """
    Calculates the tight upper bound on the treewidth of the union of two graphs.
    """
    # Example values for t_H, t_G, and k.
    # The user can change these values to their specific problem.
    t_H = 5  # treewidth of graph H
    t_G = 6  # treewidth of graph G
    k = 4    # number of shared vertices

    # The tight upper bound on the treewidth of the combined graph F is
    # max(t_H, t_G, k-1).
    if k == 0:
        # If there are no shared vertices, the treewidth is the max of the two.
        result = max(t_H, t_G)
        print(f"Since k=0, the graphs are disjoint.")
        print(f"The treewidth is max(t_H, t_G) = max({t_H}, {t_G}) = {result}")
    else:
        # For k > 0
        k_minus_1 = k - 1
        result = max(t_H, t_G, k_minus_1)
    
        # Output the variables and the final equation
        print(f"Given values:")
        print(f"t_H = {t_H}")
        print(f"t_G = {t_G}")
        print(f"k = {k}")
        print("\nThe tight upper bound for the treewidth of F is calculated as:")
        print(f"max(t_H, t_G, k-1) = max({t_H}, {t_G}, {k}-1) = max({t_H}, {t_G}, {k_minus_1}) = {result}")

solve_treewidth_bound()