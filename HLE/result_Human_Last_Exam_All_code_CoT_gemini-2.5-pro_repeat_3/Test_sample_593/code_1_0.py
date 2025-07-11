def calculate_treewidth_bound(t_H, t_G, k):
    """
    Calculates and prints the tight upper bound for the treewidth of the union of two graphs.

    Args:
        t_H (int): The treewidth of graph H.
        t_G (int): The treewidth of graph G.
        k (int): The number of shared vertices between H and G.
    """
    print(f"Given values are:")
    print(f"Treewidth of H, t_H = {t_H}")
    print(f"Treewidth of G, t_G = {t_G}")
    print(f"Number of shared vertices, k = {k}\n")

    print("The tight upper bound on the treewidth of the union graph F is given by the formula:")
    print("t_F <= max(t_H, t_G) + k\n")

    print("Substituting the given values into the formula:")
    # Using f-string to show the equation with numbers
    print(f"t_F <= max({t_H}, {t_G}) + {k}")

    max_t = max(t_H, t_G)
    print(f"First, we evaluate the maximum of the two treewidths:")
    print(f"max({t_H}, {t_G}) = {max_t}\n")

    print("Now, substitute this result back into the inequality:")
    # Using f-string to show the next step
    print(f"t_F <= {max_t} + {k}")

    upper_bound = max_t + k
    print(f"Finally, we compute the sum to find the upper bound:")
    # Using f-string to show the final result
    print(f"t_F <= {upper_bound}\n")
    print(f"So, the tight upper bound on the treewidth of F is {upper_bound}.")

# Example usage with some sample values.
# You can change these values to see the calculation for different graphs.
t_H_example = 5
t_G_example = 7
k_example = 3

calculate_treewidth_bound(t_H_example, t_G_example, k_example)