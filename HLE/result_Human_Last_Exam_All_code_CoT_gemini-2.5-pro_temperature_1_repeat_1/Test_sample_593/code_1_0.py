def solve_treewidth_union():
    """
    Calculates the tight upper bound for the treewidth of the union of two graphs.

    The problem defines a graph F as the union of two graphs H and G,
    which share k vertices.
    - t_H: treewidth of H
    - t_G: treewidth of G
    - k: number of shared vertices between H and G

    The tight upper bound on the treewidth of F (t_F) is given by the formula:
    t_F <= max(t_H, t_G, k - 1)
    """

    # Example values for the parameters
    # You can change these to test other scenarios
    t_H = 7
    t_G = 5
    k = 10

    print(f"Given parameters:")
    print(f"Treewidth of H (t_H): {t_H}")
    print(f"Treewidth of G (t_G): {t_G}")
    print(f"Number of shared vertices (k): {k}\n")

    # Calculate the components of the max function
    k_minus_1 = k - 1

    # Calculate the tight upper bound
    tw_F_bound = max(t_H, t_G, k_minus_1)

    # Output the final equation with all numbers
    print("The tight upper bound on the treewidth of F is calculated as follows:")
    print(f"t_F <= max(t_H, t_G, k - 1)")
    print(f"t_F <= max({t_H}, {t_G}, {k} - 1)")
    print(f"t_F <= max({t_H}, {t_G}, {k_minus_1})")
    print(f"t_F <= {tw_F_bound}\n")

    print(f"Thus, the tight upper bound on the treewidth of F is {tw_F_bound}.")

if __name__ == '__main__':
    solve_treewidth_union()
