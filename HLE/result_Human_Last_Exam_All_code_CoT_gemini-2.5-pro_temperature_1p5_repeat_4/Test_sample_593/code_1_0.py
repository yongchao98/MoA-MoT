def solve_treewidth_bound():
    """
    Calculates and prints the tight upper bound for the treewidth of the union of two graphs.
    """

    # We use the parameters from our tightness proof example.
    # t_H: treewidth of graph H
    # t_G: treewidth of graph G
    # k: number of vertices in the intersection of H and G
    t_H = 2
    t_G = 2
    k = 5

    # The formula for the tight upper bound on the treewidth of the union graph F is:
    # max(t_H, t_G, k-1)

    # First, calculate the value of the term k-1
    k_minus_1 = k - 1

    # Now, find the maximum of the three values to get the bound
    bound = max(t_H, t_G, k_minus_1)

    # Print the explanation and the final equation with the numbers plugged in
    print("The tight upper bound on the treewidth of F, denoted tw(F), is given by the formula:")
    print("tw(F) <= max(t_H, t_G, k-1)")
    print("\nFor the example case with t_H = {}, t_G = {}, and k = {}:".format(t_H, t_G, k))
    
    # This line prints each number in the final equation as requested.
    print(f"max({t_H}, {t_G}, {k} - 1) = max({t_H}, {t_G}, {k_minus_1}) = {bound}")

solve_treewidth_bound()