def solve_hss_multiplication_access():
    """
    Calculates the number of submatrices accessed during the multiplication
    of two matrices represented by an HSS tree of depth 4.

    The calculation is based on a recurrence relation derived from the
    block matrix multiplication algorithm applied hierarchically.

    Let T(d) be the number of submatrix accesses for a tree of depth d.
    The recurrence relation is: T(d) = 2 * T(d-1) + 12
    The base case (a leaf node, d=0) is: T(0) = 2
    """
    depth = 4
    
    # Base case: T(0) for a leaf node (depth 0)
    t_prev = 2
    print("For a tree of depth d=0 (a single leaf block):")
    print(f"T(0) = {t_prev}")
    print("-" * 20)

    # Iteratively solve the recurrence for d = 1, 2, 3, 4
    for d in range(1, depth + 1):
        # Apply the recurrence relation: T(d) = 2 * T(d-1) + 12
        t_current = 2 * t_prev + 12
        
        print(f"Calculating for a tree of depth d={d}:")
        print(f"T({d}) = 2 * T({d-1}) + 12")
        print(f"T({d}) = 2 * {t_prev} + 12 = {2 * t_prev} + 12 = {t_current}")
        print("-" * 20)
        
        t_prev = t_current

    print(f"\nFor a Hierarchical Semi-separable tree with depth {depth}, the total number of submatrices accessed is {t_prev}.")

solve_hss_multiplication_access()