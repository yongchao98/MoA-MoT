def calculate_treewidth_bound(t_H, t_G, k):
    """
    Calculates and prints the tight upper bound for the treewidth of the union of two graphs.

    Args:
        t_H (int): The treewidth of graph H.
        t_G (int): The treewidth of graph G.
        k (int): The number of vertices in the intersection of H and G.
    """
    if not all(isinstance(i, int) and i >= 0 for i in [t_H, t_G, k]):
        print("Error: Treewidths and k must be non-negative integers.")
        return

    # Determine the maximum of the two treewidths
    max_t = max(t_H, t_G)

    # Calculate the upper bound using the formula
    bound = max_t + k

    # Print the explanation and the result
    print(f"Given parameters:")
    print(f"  - Treewidth of H (t_H): {t_H}")
    print(f"  - Treewidth of G (t_G): {t_G}")
    print(f"  - Size of intersection (k): {k}")
    print("\nThe tight upper bound for the treewidth of the combined graph F is t_F <= max(t_H, t_G) + k.")
    print("\nCalculation steps:")
    # The final code outputs each number in the final equation as requested
    print(f"t_F <= max({t_H}, {t_G}) + {k}")
    print(f"t_F <= {max_t} + {k}")
    print(f"t_F <= {bound}")


# Example usage with some values for t_H, t_G, and k.
calculate_treewidth_bound(t_H=5, t_G=7, k=3)