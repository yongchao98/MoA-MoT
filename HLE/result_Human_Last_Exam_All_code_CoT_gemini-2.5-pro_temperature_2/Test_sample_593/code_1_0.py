def calculate_treewidth_bound(t_H, t_G, k):
    """
    Calculates the tight upper bound on the treewidth of a graph F,
    formed by the union of graphs H and G sharing k vertices.

    Args:
        t_H (int): The treewidth of graph H.
        t_G (int): The treewidth of graph G.
        k (int): The number of shared vertices between H and G.
    """
    if t_H < 0 or t_G < 0 or k < 0:
        print("Error: Treewidth and the number of vertices cannot be negative.")
        return

    print("The tight upper bound on the treewidth of F is given by the formula: t_F <= max(t_H, t_G, k - 1)")
    print("-" * 20)
    print("Given values:")
    print(f"t_H = {t_H}")
    print(f"t_G = {t_G}")
    print(f"k = {k}")
    print("-" * 20)

    # Calculate k-1, ensuring k>0 for k-1 to be non-negative, though the max will handle it.
    k_minus_1 = max(0, k - 1)

    # Calculate the bound
    bound = max(t_H, t_G, k_minus_1)

    print("Calculating the bound:")
    # Print the formula with the numbers plugged in
    print(f"t_F <= max({t_H}, {t_G}, {k} - 1)")
    if k > 0:
        print(f"t_F <= max({t_H}, {t_G}, {k_minus_1})")
    
    print("-" * 20)
    print(f"The tight upper bound for t_F is: {bound}")


if __name__ == '__main__':
    # Example values for t_H, t_G, and k
    example_t_H = 8
    example_t_G = 12
    example_k = 10

    calculate_treewidth_bound(example_t_H, example_t_G, example_k)
