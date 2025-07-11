import collections

def solve_knot_problem():
    """
    Calculates the specified sum for the given trefoil knot grid diagram.
    """
    n = 3
    crossings = {(0, 0), (1, 1), (2, 2)}
    lattice_points = [(i, j) for i in range(n + 1) for j in range(n + 1)]

    # Step 1: Determine the k-value for each lattice point.
    # k is the number of adjacent cells with a crossing.
    k_values = {}
    for i, j in lattice_points:
        k = 0
        # A lattice point (i,j) is a corner to cells (i-1,j-1), (i-1,j), (i,j-1), (i,j)
        if (i - 1, j - 1) in crossings: k += 1
        if (i - 1, j) in crossings: k += 1
        if (i, j - 1) in crossings: k += 1
        if (i, j) in crossings: k += 1
        k_values[(i, j)] = k

    # Step 2: Determine the winding number w(i,j) for each lattice point.
    # This is based on the known Seifert circles for this diagram.
    # Circle 1 (C1) is the boundary of cell (1,0).
    # Circle 2 (C2) is the boundary of the union of cells in rows 1 and 2.
    # w=1 if a point is on or inside either circle region, 0 otherwise.
    w_values = {}
    for i, j in lattice_points:
        w = 0
        # Check if inside or on the boundary of the region of C1 (cell (1,0))
        on_c1 = (i in [1, 2] and j in [0, 1])
        # Check if inside or on the boundary of the region of C2 (rows 1 and 2)
        on_c2 = (j >= 1)
        if on_c1 or on_c2:
            w = 1
        w_values[(i, j)] = w

    # Step 3: Compute the sum.
    # Group points by their k-value.
    points_by_k = collections.defaultdict(list)
    for point in lattice_points:
        k = k_values[point]
        if k > 0:
            points_by_k[k].append(point)

    total_sum = 0
    sum_expressions = []

    for k in sorted(points_by_k.keys()):
        # Get the winding numbers for all points with this k-value
        winds = [w_values[p] for p in points_by_k[k]]
        sum_of_winds = sum(winds)
        term = k * sum_of_winds
        total_sum += term
        
        # Format the expression for this k value
        wind_sum_str = " + ".join(map(str, winds))
        sum_expressions.append(f"{k} * ({wind_sum_str})")

    # Handle k-values for which there were no points
    for k in range(1, 5):
        if k not in points_by_k:
            sum_expressions.append(f"{k} * (0)")
            
    final_equation_str = " + ".join(sorted(sum_expressions, key=lambda x: x[0]))
    print(f"The calculation is:")
    print(f"{final_equation_str} = {total_sum}")

solve_knot_problem()