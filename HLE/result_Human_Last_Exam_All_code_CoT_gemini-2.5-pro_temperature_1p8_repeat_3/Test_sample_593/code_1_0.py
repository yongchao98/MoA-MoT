def solve_treewidth_bound():
    """
    Calculates the tight upper bound on the treewidth of a graph F
    formed by the union of two graphs H and G intersecting on k vertices.
    """
    # Placeholder values for the symbolic variables
    # t_H: treewidth of graph H
    # t_G: treewidth of graph G
    # k: number of vertices in the intersection of H and G
    t_H = 5
    t_G = 7
    k = 4

    print(f"Given symbolic values:")
    print(f"t_H (treewidth of H) = {t_H}")
    print(f"t_G (treewidth of G) = {t_G}")
    print(f"k (size of intersection) = {k}\n")
    
    # The tight upper bound on the treewidth of F is max(t_H, t_G) + k
    result = max(t_H, t_G) + k
    
    # The final equation and its result
    print("The tight upper bound is given by the formula: tw(F) <= max(t_H, t_G) + k")
    print("\nCalculating the result for the given values:")
    print(f"max({t_H}, {t_G}) + {k} = {result}")

solve_treewidth_bound()