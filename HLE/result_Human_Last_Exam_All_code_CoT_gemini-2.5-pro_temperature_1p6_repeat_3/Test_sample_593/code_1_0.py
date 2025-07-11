def solve():
    """
    Calculates the tight upper bound on the treewidth of the union of two graphs.

    The problem states:
    Let H be a graph of treewidth t_H, and G be a graph of treewidth t_G.
    Let V(H) intersect V(G) on a set U of k vertices.
    F is the union of H and G.

    The tight upper bound for the treewidth of F is max(t_H, t_G) + k.
    """
    # Example values for t_H, t_G, and k
    t_H = 7
    t_G = 5
    k = 4

    # The formula for the tight upper bound on the treewidth of F
    bound = max(t_H, t_G) + k
    
    # Output the explanation and the calculation
    print("Let t_H be the treewidth of graph H, t_G be the treewidth of graph G, and k be the number of shared vertices.")
    print(f"Given values: t_H = {t_H}, t_G = {t_G}, k = {k}")
    print("\nThe formula for the tight upper bound on the treewidth of the combined graph F = H U G is:")
    print("tw(F) <= max(t_H, t_G) + k")
    print("\nCalculation:")
    print(f"max({t_H}, {t_G}) + {k} = {max(t_H, t_G)} + {k} = {bound}")
    print(f"\nThe tight upper bound is: {bound}")

solve()