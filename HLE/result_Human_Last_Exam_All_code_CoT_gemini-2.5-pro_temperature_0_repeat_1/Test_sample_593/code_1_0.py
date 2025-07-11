def calculate_treewidth_bound(t_H, t_G, k):
    """
    Calculates the tight upper bound on the treewidth of a graph F
    formed by the union of two graphs H and G.

    Args:
        t_H (int): The treewidth of graph H.
        t_G (int): The treewidth of graph G.
        k (int): The number of shared vertices between H and G.
    """
    if not all(isinstance(i, int) and i >= 0 for i in [t_H, t_G, k]):
        print("Error: Treewidths and k must be non-negative integers.")
        return

    # The tight upper bound is max(t_H, t_G, k-1)
    # We handle the case k=0, where k-1 would be -1. Treewidth is non-negative.
    k_term = max(0, k - 1)
    
    bound = max(t_H, t_G, k_term)

    # Output the explanation and the result
    print(f"Given values:")
    print(f"  Treewidth of H (t_H): {t_H}")
    print(f"  Treewidth of G (t_G): {t_G}")
    print(f"  Number of shared vertices (k): {k}")
    print("\nThe tight upper bound on the treewidth of the combined graph F is given by the formula:")
    print("  tw(F) <= max(t_H, t_G, k-1)")
    print("\nCalculation:")
    print(f"  max({t_H}, {t_G}, {k}-1) = max({t_H}, {t_G}, {k-1}) = {bound}")
    print(f"\nThe calculated tight upper bound is: {bound}")

# Example usage with some sample values.
# You can change these values to see the result for different scenarios.
example_t_H = 5
example_t_G = 7
example_k = 4

calculate_treewidth_bound(example_t_H, example_t_G, example_k)