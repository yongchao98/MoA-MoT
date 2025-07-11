def solve_treewidth_bound():
    """
    Calculates and prints the tight upper bound on the treewidth of a graph F
    formed by the union of two graphs H and G sharing k vertices.
    """

    # Example values for the treewidths of H and G, and the number of common vertices.
    # You can change these values to see the result for a different case.
    t_H = 8
    t_G = 10
    k = 12

    print(f"Given values:")
    print(f"Treewidth of H, t_H = {t_H}")
    print(f"Treewidth of G, t_G = {t_G}")
    print(f"Number of common vertices, k = {k}\n")

    # The formula for the tight upper bound on the treewidth of F is max(t_H, t_G, k-1).

    # Calculate the components of the max function
    val1 = t_H
    val2 = t_G
    val3 = k - 1

    # Calculate the final bound using the max function
    treewidth_F_bound = max(val1, val2, val3)

    # Print the result showing the full equation as requested
    print(f"The tight upper bound for the treewidth of F is given by the formula: max(t_H, t_G, k-1)")
    print(f"Calculation:")
    print(f"Bound = max({val1}, {val2}, {k}-1)")
    print(f"Bound = max({val1}, {val2}, {val3}) = {treewidth_F_bound}")

# Execute the function to see the output
solve_treewidth_bound()