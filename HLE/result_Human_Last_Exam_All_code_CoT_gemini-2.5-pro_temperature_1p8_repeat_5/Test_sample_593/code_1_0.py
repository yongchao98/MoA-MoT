# User-defined parameters for the treewidth problem.
# You can change these values to explore different scenarios.

# t_H: Treewidth of graph H
t_H = 4

# t_G: Treewidth of graph G
t_G = 5

# k: Number of vertices in the intersection of H and G
k = 8


def solve_treewidth_union_bound(t_H, t_G, k):
    """
    Calculates and prints the tight upper bound on the treewidth of the union of two graphs.

    The tight upper bound for the treewidth of graph F = H U G is given by the formula:
    tw(F) <= max(t_H, t_G, k-1)

    Args:
        t_H (int): The treewidth of graph H.
        t_G (int): The treewidth of graph G.
        k (int): The number of vertices in the intersection V(H) intersect V(G).
    """

    # The problem assumes k >= 1 since an empty intersection is trivial.
    # The 'no isolated vertices' constraint implies t_H >= 1 and t_G >= 1 for k > 1.
    # However, the formula holds even without these constraints.
    if k <= 0:
        # If the intersection is empty, the treewidth is the max of the individual treewidths.
        # Our formula max(t_H, t_G, -1) still works if max is defined for this case.
        # But to be clear, we handle k=0 as a separate case.
        bound = max(t_H, t_G)
        print("The graphs are disjoint (k=0).")
        print(f"The treewidth of the union is max(t_H, t_G) = max({t_H}, {t_G}) = {bound}")
        return

    # Calculate k-1 term
    k_minus_1 = k - 1
    
    # Calculate the final bound
    bound = max(t_H, t_G, k_minus_1)

    # Print the explanation and the result
    print("Let t_H be the treewidth of graph H.")
    print("Let t_G be the treewidth of graph G.")
    print("Let k be the number of vertices in the intersection of H and G.")
    
    print("\nA tight upper bound for the treewidth of the union graph F = H U G is max(t_H, t_G, k-1).\n")

    print(f"For the given values:")
    print(f"t_H = {t_H}")
    print(f"t_G = {t_G}")
    print(f"k = {k}\n")
    
    print("The final equation is:")
    # Showing each number in the final equation as requested
    print(f"max({t_H}, {t_G}, {k} - 1) = max({t_H}, {t_G}, {k_minus_1}) = {bound}")

solve_treewidth_union_bound(t_H, t_G, k)