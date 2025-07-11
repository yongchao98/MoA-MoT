def calculate_treewidth_bound():
    """
    Calculates and prints a tight upper bound on the treewidth of the union of two graphs.
    """
    #
    # Define the treewidths of graphs H and G, and the number of shared vertices k.
    # These are example values. You can change them to any non-negative integer.
    #
    # t_H: Treewidth of graph H
    t_H = 5
    # t_G: Treewidth of graph G
    t_G = 3
    # k: Number of shared vertices in V(H) intersect V(G)
    k = 7

    # The problem states neither H nor G has isolated vertices.
    # We must have t_H >= 0, t_G >= 0, and k >= 0.

    # The formula for the tight upper bound on the treewidth of F = H U G is max(t_H, t_G, k-1).
    # We calculate k-1 first. Note that if k=0, this is -1, which is fine as max will handle it.
    k_minus_1 = k - 1

    # The treewidth of F is the maximum of the three values.
    treewidth_F = max(t_H, t_G, k_minus_1)

    print(f"Given the treewidth of H as t_H = {t_H}, the treewidth of G as t_G = {t_G}, and the number of shared vertices as k = {k}:")
    print("A tight upper bound for the treewidth of the combined graph F is calculated using the formula:")
    print("tw(F) <= max(t_H, t_G, k-1)")
    print()
    print("Calculating the values:")
    print(f"t_H = {t_H}")
    print(f"t_G = {t_G}")
    print(f"k-1 = {k} - 1 = {k_minus_1}")
    print()
    print("The final bound is:")
    print(f"tw(F) <= max({t_H}, {t_G}, {k_minus_1}) = {treewidth_F}")

# Execute the function to see the result.
calculate_treewidth_bound()