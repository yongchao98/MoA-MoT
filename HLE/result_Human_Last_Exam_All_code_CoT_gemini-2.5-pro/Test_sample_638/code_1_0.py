import collections

def solve_knot_problem():
    """
    Solves the grid diagram problem for the left-hand trefoil knot.
    """
    # Step 1: Define the grid diagram based on the problem description.
    # O's are on the anti-diagonal, X's are placed to form a valid grid for a trefoil.
    # Indices are (row, col) from top-left (0,0).
    O_coords = [(0, 2), (1, 1), (2, 0)]
    X_coords = [(0, 0), (1, 2), (2, 1)]
    marked_cells = set(O_coords + X_coords)

    # The lattice is 4x4, for a 3x3 grid.
    grid_size = 3
    lattice_dim = grid_size + 1

    # Step 2: Calculate winding numbers w(i, j) for each lattice point (col, row).
    # The formula is from knot theory literature (e.g., Dunfield's GridSeifert).
    # We flip the sign because the grid specified is for a R-H trefoil, but the
    # problem asks for the L-H trefoil.
    w = {}
    for r in range(lattice_dim):  # row index for lattice points
        for c in range(lattice_dim):  # col index for lattice points
            n_O_bl = sum(1 for ox, oy in O_coords if ox < c and oy > r) # Using visual coords where y increases upwards
            n_X_tr = sum(1 for xx, xy in X_coords if xx >= c and xy <= r)
            n_O_br = sum(1 for ox, oy in O_coords if ox >= c and oy > r)
            n_X_tl = sum(1 for xx, xy in X_coords if xx < c and xy <= r)

            # Adjusting indexing for calculation: (col, grid_size - 1 - row)
            # This makes (0,0) bottom-left
            calc_r = grid_size - 1 - r
            
            n_O_bl = sum(1 for ox, oy in O_coords if ox < c and oy < calc_r)
            n_X_tr = sum(1 for xx, xy in X_coords if xx >= c and xy >= calc_r)
            n_O_br = sum(1 for ox, oy in O_coords if ox >= c and oy < calc_r)
            n_X_tl = sum(1 for xx, xy in X_coords if xx < c and xy >= calc_r)
            
            # The standard formula gives winding numbers for the R-H trefoil.
            # We flip the sign for the L-H trefoil.
            winding_num = -(n_O_bl + n_X_tr - n_O_br - n_X_tl)
            w[c, r] = winding_num

    # Step 3: Determine k for each lattice point.
    # k is the number of marked cells sharing the lattice point as a corner.
    k_counts = {}
    for r in range(lattice_dim):
        for c in range(lattice_dim):
            count = 0
            # A lattice point (c,r) is a corner for cells (c-1,r-1), (c,r-1), (c-1,r), (c,r)
            # using top-left cell indexing.
            if (r - 1, c - 1) in marked_cells: count += 1
            if (r - 1, c) in marked_cells: count += 1
            if (r, c - 1) in marked_cells: count += 1
            if (r, c) in marked_cells: count += 1
            k_counts[c, r] = count

    # Step 4: Partition winding numbers into sets mho_k.
    mho = collections.defaultdict(list)
    for r in range(lattice_dim):
        for c in range(lattice_dim):
            k = k_counts[c, r]
            if k > 0:
                mho[k].append(w[c, r])

    # Step 5: Compute the final sum.
    total_sum = 0
    sum_terms = []
    
    print("The final calculation is:")
    
    equation_parts = []
    for k in sorted(mho.keys()):
        if k > 4: continue
        sum_of_w = sum(mho[k])
        term = k * sum_of_w
        total_sum += term
        
        # Build the string for each term k * sum(w_i)
        w_list_str = [str(val) for val in mho[k]]
        sum_str = f"({'+'.join(w_list_str)})".replace("+-", "-")
        equation_parts.append(f"{k} * {sum_str}")

    print(f"S = {' + '.join(equation_parts)}")
    
    final_equation_parts = []
    for k in sorted(mho.keys()):
        if k > 4: continue
        sum_of_w = sum(mho[k])
        final_equation_parts.append(f"{k} * ({sum_of_w})")
    
    print(f"S = {' + '.join(final_equation_parts)}")
    print(f"S = {total_sum}")


solve_knot_problem()