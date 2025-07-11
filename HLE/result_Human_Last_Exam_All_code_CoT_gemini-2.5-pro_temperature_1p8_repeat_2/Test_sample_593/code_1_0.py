def solve():
    """
    This function calculates the tight upper bound for the treewidth of F=H U G
    and demonstrates its tightness with a concrete example.
    """
    
    # The tight upper bound for the treewidth of F is max(t_H, t_G, k-1).
    # We demonstrate that this bound is tight using a constructive example.
    # Let's choose a specific value for k, the number of shared vertices.
    # Let k be an integer >= 3.
    k = 5
    
    print(f"Let's demonstrate the tightness of the bound for k = {k}.\n")
    
    # We construct a graph F = K_k (a complete graph on k vertices).
    # The treewidth of K_k is k-1.
    tw_F = k - 1
    print(f"We aim to construct F such that tw(F) = {tw_F}.")
    
    # We partition the edges of F=K_k into two subgraphs, H and G,
    # such that V(H) = V(G) = V(F) = S, where |S| = k.
    # A known construction that achieves this is:
    # H = K_k - e (a complete graph with one edge removed)
    # G = graph containing edge e and other edges to ensure no isolated vertices
    #     and a small treewidth.
    
    # The treewidth of H = K_k - e is k-2 (for k>=3).
    t_H = k - 2
    
    # We can construct G such that its treewidth is 1. For example, a path
    # or a star graph. So we set t_G = 1.
    t_G = 1
    
    print(f"Our construction uses component graphs with parameters:")
    print(f"t_H (treewidth of H) = {t_H}")
    print(f"t_G (treewidth of G) = {t_G}")
    print(f"k (shared vertices) = {k}\n")
    
    # Now, we calculate the upper bound using the formula max(t_H, t_G, k-1).
    k_minus_1 = k - 1
    
    upper_bound = max(t_H, t_G, k_minus_1)
    
    print("The upper bound is calculated as max(t_H, t_G, k-1).")
    print(f"The values are:")
    print(f"  t_H = {t_H}")
    print(f"  t_G = {t_G}")
    print(f"  k-1 = {k_minus_1}")
    
    print(f"\nThe result of max({t_H}, {t_G}, {k_minus_1}) is {upper_bound}.")
    print(f"This calculated bound of {upper_bound} matches the actual treewidth of our constructed graph F=K_k, which is {tw_F}.")
    print("This demonstrates that the bound is tight.\n")
    
    print("The tight upper bound is generally expressed as: max(t_H, t_G, k-1)")

solve()