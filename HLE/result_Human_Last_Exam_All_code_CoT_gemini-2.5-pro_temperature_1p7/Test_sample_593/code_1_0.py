def solve_treewidth_bound():
    """
    This function calculates and prints the tight upper bound on the treewidth of
    a graph F formed by the union of two graphs H and G.

    The user should define the following variables:
    t_H: The treewidth of graph H.
    t_G: The treewidth of graph G.
    k: The number of vertices in the intersection of H and G.
    """
    # Example values.
    # You can change these to reflect a specific problem instance.
    t_H = 4
    t_G = 5
    k = 8

    # The tight upper bound is max(t_H, t_G, k-1)
    bound = max(t_H, t_G, k - 1)

    print("The tight upper bound on the treewidth of F is given by the formula: max(t_H, t_G, k - 1).")
    print(f"For the given values:")
    print(f"t_H = {t_H}")
    print(f"t_G = {t_G}")
    print(f"k = {k}")
    print("\nThe calculation is:")
    print(f"max({t_H}, {t_G}, {k} - 1) = max({t_H}, {t_G}, {k - 1}) = {bound}")
    print(f"\nThe treewidth of F is at most {bound}.")

solve_treewidth_bound()