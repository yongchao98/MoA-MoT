def solve():
    """
    Calculates the tight upper bound on the treewidth of the union of two graphs.
    """
    # Let's use some example values for t_H, t_G, and k.
    # Treewidth of graph H
    t_H = 10
    # Treewidth of graph G
    t_G = 8
    # Number of vertices in the intersection of H and G
    k = 4

    # The formula for the tight upper bound on the treewidth of F is max(t_H, t_G) + k.
    max_treewidth = max(t_H, t_G)
    upper_bound = max_treewidth + k
    
    # Print the equation step by step.
    print(f"Let t_H be the treewidth of H, t_G be the treewidth of G, and k be the number of common vertices.")
    print(f"A tight upper bound for the treewidth of F = H U G is given by the formula: max(t_H, t_G) + k")
    print(f"Given values: t_H = {t_H}, t_G = {t_G}, k = {k}")
    print(f"The calculation is: max({t_H}, {t_G}) + {k} = {max_treewidth} + {k} = {upper_bound}")
    print(f"So, the tight upper bound on the treewidth of F is {upper_bound}.")

solve()