def calculate_treewidth_upper_bound(t_H, t_G, k):
    """
    Calculates the tight upper bound on the treewidth of the union of two graphs.

    Args:
        t_H (int): The treewidth of graph H.
        t_G (int): The treewidth of graph G.
        k (int): The number of vertices in the intersection of H and G.
    """
    if not all(isinstance(i, int) and i >= 0 for i in [t_H, t_G, k]):
        print("Error: Treewidths and k must be non-negative integers.")
        return

    # The tight upper bound is max(t_H, t_G, k-1)
    # Ensure k-1 is at least 0 if k is 0, though k >= 1 is more practical.
    k_minus_1 = max(0, k - 1)
    upper_bound = max(t_H, t_G, k_minus_1)

    print("A tight upper bound for the treewidth of the combined graph F is given by the formula: max(t_H, t_G, k-1)")
    print("\nGiven the input values:")
    print(f"  Treewidth of H (t_H): {t_H}")
    print(f"  Treewidth of G (t_G): {t_G}")
    print(f"  Number of shared vertices (k): {k}")
    
    print(f"\nThe calculation is:")
    print(f"  max({t_H}, {t_G}, {k}-1) = max({t_H}, {t_G}, {k_minus_1})")
    print(f"  Result: {upper_bound}")

# --- Example Usage ---
# You can change these values to test with different graphs.
example_t_H = 5
example_t_G = 7
example_k = 6

calculate_treewidth_upper_bound(example_t_H, example_t_G, example_k)