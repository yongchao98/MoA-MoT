def calculate_treewidth_bound(t_H, t_G, k):
    """
    Calculates the tight upper bound on the treewidth of F = H U G.

    Args:
        t_H (int): The treewidth of graph H.
        t_G (int): The treewidth of graph G.
        k (int): The number of vertices in the intersection of H and G.
    """
    # The problem statement implies H and G are not trivial graphs,
    # so their treewidth must be at least 1, and k must be at least 1
    # if the intersection is non-empty. We can add checks for valid inputs.
    if t_H < 0 or t_G < 0 or k < 0:
        print("Error: Treewidth and k cannot be negative.")
        return

    # The formula for the tight upper bound is max(t_H, t_G, k-1)
    k_minus_1 = k - 1
    bound = max(t_H, t_G, k_minus_1)

    # Print the explanation and the final equation with numbers
    print("A tight upper bound for the treewidth of the union graph F is given by the formula:")
    print("tw(F) <= max(t_H, t_G, k - 1)")
    print("\nFor the given values:")
    print(f"t_H = {t_H}")
    print(f"t_G = {t_G}")
    print(f"k = {k}")
    print("\nThe calculation is:")
    print(f"max({t_H}, {t_G}, {k} - 1) = max({t_H}, {t_G}, {k_minus_1}) = {bound}")
    print(f"\nThe tight upper bound is: {bound}")

# Example Usage:
# Let's assume H has treewidth 5, G has treewidth 4,
# and they intersect on 7 vertices.
example_t_H = 5
example_t_G = 4
example_k = 7
calculate_treewidth_bound(example_t_H, example_t_G, example_k)