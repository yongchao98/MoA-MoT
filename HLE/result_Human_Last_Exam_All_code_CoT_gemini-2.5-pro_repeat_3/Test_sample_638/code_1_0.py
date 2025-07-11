import collections

def solve_knot_problem():
    """
    Solves the grid diagram winding number problem for the left-hand trefoil.
    """
    n = 3

    # Step 1: Define the grid diagram for the left-hand trefoil.
    # O's are on the diagonal (p_O = id).
    # X's correspond to the permutation p_X = (1 3 2) for the left-hand trefoil.
    # We use 1-based indexing for rows and columns, as is standard.
    o_locs = {(1, 1), (2, 2), (3, 3)}
    p_x = {1: 3, 2: 1, 3: 2}
    x_locs = set(p_x.items())
    markers = o_locs.union(x_locs)

    # Step 2: Calculate winding numbers w(i,j) for each lattice point (i,j).
    # Lattice points (i,j) are 0-indexed, from 0 to n.
    w_matrix = [[0] * (n + 1) for _ in range(n + 1)]
    for i in range(n + 1):
        for j in range(n + 1):
            # Formula: w(i,j) = |{k | k>i, p_X(k)<j}| - |{k | k<i, p_X(k)>j}|
            # Here, i and j are 0-based lattice coordinates.
            # k and p_x(k) are 1-based column/row indices.
            
            term1_count = 0
            for k in range(1, n + 1):
                if k > i and p_x[k] < j:
                    term1_count += 1
            
            term2_count = 0
            for k in range(1, n + 1):
                if k < i and p_x[k] > j:
                    term2_count += 1
            
            w_matrix[j][i] = term1_count - term2_count

    # Step 3: Determine the sets mho_k.
    mho = collections.defaultdict(list)
    for i in range(n + 1):
        for j in range(n + 1):
            # Count markers in the four cells meeting at lattice point (i,j).
            # The four cells are (i,j), (i+1,j), (i,j+1), (i+1,j+1) in 1-based coordinates.
            k_count = 0
            # Cells are indexed by (column, row)
            adjacent_cells = [(i, j), (i + 1, j), (i, j + 1), (i + 1, j + 1)]
            for c, r in adjacent_cells:
                if (c, r) in markers:
                    k_count += 1
            
            if k_count > 0:
                mho[k_count].append(w_matrix[j][i])

    # Step 4: Compute and print the final sum.
    total_sum = 0
    equation_parts = []
    for k in range(1, 5):
        if not mho[k]:
            sum_w = 0
            w_str = "(0)"
        else:
            sum_w = sum(mho[k])
            w_str = "(" + " + ".join(map(str, mho[k])).replace("+ -", "- ") + ")"
        
        term = k * sum_w
        total_sum += term
        equation_parts.append(f"{k} * {w_str}")

    print("The final calculation is:")
    print(" + ".join(equation_parts) + f" = {total_sum}")

solve_knot_problem()