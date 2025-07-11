def calculate_treewidth_union(t_H, t_G, k):
    """
    Calculates and displays the tight upper bound for the treewidth of the union of two graphs.

    The problem describes a graph F formed by the union of two graphs H and G,
    which intersect on a set of k vertices. A key result in graph theory states
    that for such a composition, the treewidth of the resulting graph F is the
    maximum of the treewidths of the original graphs H and G.

    tw(F) = max(tw(H), tw(G))

    This bound is tight and, interestingly, does not depend on the number of
    shared vertices, k.

    Args:
      t_H (int): The treewidth of graph H.
      t_G (int): The treewidth of graph G.
      k (int): The number of vertices in the intersection of H and G.
    """
    
    # The tight upper bound is the maximum of the two treewidths.
    bound = max(t_H, t_G)

    print("The tight upper bound on the treewidth of F, tw(F), is given by the formula:")
    print("tw(F) = max(t_H, t_G)")
    print("\nFor the given values:")
    print(f"t_H = {t_H}")
    print(f"t_G = {t_G}")
    # The variable k is not used in the final formula, but we acknowledge it.
    print(f"k = {k}")
    
    print("\nThe final equation is:")
    # The instruction requires printing each number in the final equation.
    print(f"max({t_H}, {t_G}) = {bound}")

# Example usage with arbitrary values for t_H, t_G, and k.
# You can change these values to test with other numbers.
example_t_H = 5
example_t_G = 7
example_k = 3

calculate_treewidth_union(example_t_H, example_t_G, example_k)