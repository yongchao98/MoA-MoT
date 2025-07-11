import sys

def solve_treewidth_union():
    """
    Calculates the tight upper bound on the treewidth of the union of two graphs.

    The problem defines a graph F as the union of two graphs H and G,
    which intersect on a set of k vertices.
    - H has treewidth t_H.
    - G has treewidth t_G.
    - The number of shared vertices is k.

    The tight upper bound for the treewidth of the resulting graph F is
    given by the formula: max(t_H, t_G, k-1).
    """

    # Example values for t_H, t_G, and k.
    # You can change these to see the result for different graphs.
    t_H = 7
    t_G = 5
    k = 6

    # The treewidth of a k-clique is k-1. This is the third term in our max function.
    k_minus_1 = k - 1

    # The upper bound is the maximum of the three values.
    treewidth_F_bound = max(t_H, t_G, k_minus_1)

    # Output the explanation and the final equation with all numbers.
    print("Let t_H be the treewidth of graph H, t_G be the treewidth of graph G, and k be the number of shared vertices.")
    print(f"For the given values, t_H = {t_H}, t_G = {t_G}, and k = {k}.")
    print("\nThe tight upper bound on the treewidth of the combined graph F is given by the formula: max(t_H, t_G, k-1).")
    print("\nCalculating the result:")
    # Using f-string to format the final output equation clearly.
    print(f"max({t_H}, {t_G}, {k}-1) = max({t_H}, {t_G}, {k_minus_1}) = {treewidth_F_bound}")


if __name__ == '__main__':
    solve_treewidth_union()