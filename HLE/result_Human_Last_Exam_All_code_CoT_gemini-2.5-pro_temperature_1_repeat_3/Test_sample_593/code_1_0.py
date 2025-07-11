def solve_treewidth_bound():
    """
    Calculates and prints the tight upper bound on the treewidth of the
    union of two graphs H and G.
    """
    # Example values for the treewidths of H and G, and the size of their intersection.
    # Let t_H be the treewidth of graph H.
    t_H = 7
    # Let t_G be the treewidth of graph G.
    t_G = 10
    # Let k be the number of vertices in the intersection of H and G.
    k = 4

    # The tight upper bound for the treewidth of the combined graph F is given by the formula:
    # max(t_H, t_G) + k
    
    # Calculate the maximum of the two treewidths
    max_t = max(t_H, t_G)
    
    # Calculate the final bound
    result = max_t + k
    
    # Print the equation and the final answer, showing each number.
    print(f"Let t_H = {t_H}, t_G = {t_G}, and k = {k}")
    print("The tight upper bound on the treewidth of F is max(t_H, t_G) + k.")
    print(f"max({t_H}, {t_G}) + {k} = {max_t} + {k} = {result}")

solve_treewidth_bound()
