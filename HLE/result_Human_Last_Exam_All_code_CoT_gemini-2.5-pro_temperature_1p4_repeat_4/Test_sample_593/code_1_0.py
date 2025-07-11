def calculate_treewidth_bound(t_H, t_G, k):
    """
    Calculates the tight upper bound on the treewidth of the union of two graphs.

    Args:
        t_H (int): The treewidth of graph H.
        t_G (int): The treewidth of graph G.
        k (int): The number of vertices in the intersection of H and G.
    """
    
    # The tight upper bound is max(t_H, t_G, k-1)
    bound = max(t_H, t_G, k - 1)
    
    # We print the result based on the tightness example discussed above.
    # In our example, H = C5, G = C5 (complement in K5), F = K5
    # t_H = 2, t_G = 2, k = 5
    # The actual treewidth of F = K5 is 4.
    
    print(f"Given values from the tightness example:")
    print(f"t_H (treewidth of H = C5) = {t_H}")
    print(f"t_G (treewidth of G = complement of C5 in K5) = {t_G}")
    print(f"k (number of shared vertices) = {k}")
    print("-" * 20)
    print("The tight upper bound is calculated as: max(t_H, t_G, k-1)")
    print(f"Bound = max({t_H}, {t_G}, {k}-1) = max({t_H}, {t_G}, {k-1}) = {bound}")
    print(f"This matches the actual treewidth of the resulting graph F = K5, which is {bound}.")

# Example values from the tightness proof
t_H_example = 2
t_G_example = 2
k_example = 5

calculate_treewidth_bound(t_H_example, t_G_example, k_example)