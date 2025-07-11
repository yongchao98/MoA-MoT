def solve():
    """
    Calculates the tight upper bound on the treewidth of the union of two graphs.
    """
    # Define symbolic variables for the explanation.
    # For a concrete calculation, these would be numerical values.
    # For example: t_H = 5, t_G = 7, k = 4
    t_H = 5
    t_G = 7
    k = 4

    # The tight upper bound is given by the formula: min(t_H, t_G) + k
    upper_bound = min(t_H, t_G) + k

    # Output the explanation and the final equation with values
    print("Let t_F be the treewidth of the resulting graph F.")
    print("A tight upper bound for t_F is given by the expression: min(t_H, t_G) + k")
    print("\nFor the example values:")
    print(f"t_H = {t_H}")
    print(f"t_G = {t_G}")
    print(f"k = {k}")
    print("\nThe calculation is:")
    print(f"t_F <= min({t_H}, {t_G}) + {k} = {min(t_H, t_G)} + {k} = {upper_bound}")

solve()