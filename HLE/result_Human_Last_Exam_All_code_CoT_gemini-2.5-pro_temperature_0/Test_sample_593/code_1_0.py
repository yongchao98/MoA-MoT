def solve_treewidth_bound():
    """
    Calculates and prints the tight upper bound on the treewidth of F = H U G.
    """
    # Example values for the treewidths of H and G, and the number of shared vertices k.
    t_H = 7
    t_G = 5
    k = 4

    # The tight upper bound is max(t_H, t_G) + k.
    max_t = max(t_H, t_G)
    upper_bound = max_t + k

    print(f"Given:")
    print(f"  Treewidth of H, t_H = {t_H}")
    print(f"  Treewidth of G, t_G = {t_G}")
    print(f"  Number of shared vertices, k = {k}")
    print(f"\nThe tight upper bound for the treewidth of F is calculated as:")
    print(f"  tw(F) <= max(t_H, t_G) + k")
    print(f"  tw(F) <= max({t_H}, {t_G}) + {k}")
    print(f"  tw(F) <= {max_t} + {k}")
    print(f"  tw(F) <= {upper_bound}")

solve_treewidth_bound()