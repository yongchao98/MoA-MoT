def calculate_treewidth_bound(t_H, t_G, k):
    """
    Calculates the tight upper bound on the treewidth of the union of two graphs.

    Args:
        t_H (int): The treewidth of graph H.
        t_G (int): The treewidth of graph G.
        k (int): The number of shared vertices between H and G.
    """
    
    # The tight upper bound formula is max(t_H, t_G) + k
    bound = max(t_H, t_G) + k
    
    # Output the explanation of the final equation
    print(f"The treewidth of graph F, t_F, is bounded by:")
    print(f"t_F <= max(t_H, t_G) + k")
    print(f"Given values are t_H = {t_H}, t_G = {t_G}, and k = {k}.")
    print(f"Plugging in the values:")
    print(f"t_F <= max({t_H}, {t_G}) + {k}")
    print(f"t_F <= {max(t_H, t_G)} + {k}")
    print(f"t_F <= {bound}")
    print(f"\nThe tight upper bound on the treewidth of F is: {bound}")

if __name__ == '__main__':
    # You can change these values to test with different treewidths and intersection sizes.
    # Example values:
    t_H_val = 4
    t_G_val = 5
    k_val = 3
    
    calculate_treewidth_bound(t_H_val, t_G_val, k_val)