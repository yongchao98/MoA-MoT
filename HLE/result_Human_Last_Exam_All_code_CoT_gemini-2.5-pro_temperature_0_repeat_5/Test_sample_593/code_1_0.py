def solve_treewidth_bound():
    """
    Calculates and explains the tight upper bound on the treewidth of the union of two graphs.
    """
    # Example values for the treewidths of H and G, and the size of their intersection.
    # The user can change these values to see the result for their specific case.
    t_H = 5
    t_G = 4
    k = 3

    # The formula for the tight upper bound
    # tw(F) <= max(tw(H), tw(G)) + k
    bound = max(t_H, t_G) + k
    
    print("The problem is to find a tight upper bound on the treewidth of a graph F,")
    print("where F is the union of two graphs H and G that intersect on k vertices.")
    print("")
    print("Let t_H be the treewidth of H.")
    print("Let t_G be the treewidth of G.")
    print("Let k be the number of vertices in the intersection of H and G.")
    print("")
    print("A tight upper bound on the treewidth of F, t_F, is given by the formula:")
    print("t_F <= max(t_H, t_G) + k")
    print("")
    print("Using the example values:")
    print(f"t_H = {t_H}")
    print(f"t_G = {t_G}")
    print(f"k = {k}")
    print("")
    print("The final equation is:")
    print(f"max({t_H}, {t_G}) + {k} = {max(t_H, t_G)} + {k} = {bound}")

solve_treewidth_bound()