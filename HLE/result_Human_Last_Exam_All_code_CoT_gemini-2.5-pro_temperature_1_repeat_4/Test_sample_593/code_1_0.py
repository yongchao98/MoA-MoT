def solve_treewidth_union():
    """
    Calculates the tight upper bound on the treewidth of the union of two graphs.
    """
    # Please provide the treewidth of graph H, treewidth of graph G,
    # and the number of vertices in their intersection.
    # Example values are provided below.
    t_H = 4  # Treewidth of graph H
    t_G = 5  # Treewidth of graph G
    k = 3    # Number of vertices in the intersection V(H) intersect V(G)

    # The tight upper bound is max(t_H, t_G, k-1).
    # We must ensure k-1 is not negative, though k must be >= 0.
    # If k=0 (disjoint union), the intersection term is not applicable in the same way,
    # and the treewidth is max(t_H, t_G). The formula with k-1=-1 still works.
    if k < 0:
        print("Error: The number of intersecting vertices (k) cannot be negative.")
        return

    intersection_term = k - 1
    
    # Calculate the maximum of the three values
    upper_bound = max(t_H, t_G, intersection_term)

    # Output the result in a descriptive format
    print(f"Given values:")
    print(f"  Treewidth of H (t_H): {t_H}")
    print(f"  Treewidth of G (t_G): {t_G}")
    print(f"  Size of intersection (k): {k}")
    print("\nThe formula for the tight upper bound on the treewidth of the union F = H U G is:")
    print(f"  tw(F) <= max(t_H, t_G, k - 1)")
    print("\nCalculation:")
    print(f"  max({t_H}, {t_G}, {k} - 1) = max({t_H}, {t_G}, {intersection_term}) = {upper_bound}")
    print(f"\nThe tight upper bound is: {upper_bound}")

if __name__ == "__main__":
    solve_treewidth_union()