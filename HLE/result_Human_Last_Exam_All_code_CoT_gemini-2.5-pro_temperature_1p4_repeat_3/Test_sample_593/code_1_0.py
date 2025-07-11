def calculate_treewidth_bound(t_H, t_G, k):
    """
    Calculates and prints the tight upper bound for the treewidth of F.

    Args:
      t_H: The treewidth of graph H.
      t_G: The treewidth of graph G.
      k: The number of vertices in the intersection of H and G.
    """
    if not all(isinstance(i, int) and i >= 0 for i in [t_H, t_G, k]):
        print("Error: Treewidths and k must be non-negative integers.")
        return

    bound = max(t_H, t_G) + k

    print("The tight upper bound on the treewidth of F is given by the expression: max(t_H, t_G) + k")
    print("\nFor the specific values:")
    print(f"t_H = {t_H}")
    print(f"t_G = {t_G}")
    print(f"k = {k}")
    print("\nThe calculation is:")
    print(f"max({t_H}, {t_G}) + {k} = {bound}")

# Example usage with some placeholder values
# In a real scenario, these values would be determined by the specific graphs H and G.
example_t_H = 7
example_t_G = 5
example_k = 4

calculate_treewidth_bound(example_t_H, example_t_G, example_k)