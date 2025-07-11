def solve_treewidth_union():
    """
    Calculates and explains the tight upper bound for the treewidth of the union of two graphs.
    """
    # Let t_H be the treewidth of graph H.
    # We use an example value for demonstration.
    t_H = 7

    # Let t_G be the treewidth of graph G.
    # We use an example value for demonstration.
    t_G = 5

    # Let k be the number of vertices in the intersection of H and G.
    # We use an example value for demonstration.
    k = 4

    # The tight upper bound for the treewidth of the union graph F is max(t_H, t_G) + k.
    max_treewidth = max(t_H, t_G)
    final_bound = max_treewidth + k

    print("Let t_H be the treewidth of graph H.")
    print("Let t_G be the treewidth of graph G.")
    print("Let k be the number of vertices in the intersection V(H) and V(G).")
    print("\nA tight upper bound on the treewidth of the combined graph F (tw(F)) is given by the formula:")
    print("tw(F) <= max(t_H, t_G) + k")
    print("\nUsing the example values:")
    print(f"t_H = {t_H}")
    print(f"t_G = {t_G}")
    print(f"k = {k}")
    print("\nThe calculation for the bound is:")
    print(f"tw(F) <= max({t_H}, {t_G}) + {k}")
    print(f"       = {max_treewidth} + {k}")
    print(f"       = {final_bound}")

solve_treewidth_union()