def solve_braid_index():
    """
    Calculates the braid index for a knot defined by a grid diagram.
    """
    n = 7
    o_pos = [(1,1), (2,7), (3,4), (4,5), (5,3), (6,6), (7,2)]
    x_pos = [(1,2), (2,6), (3,3), (4,1), (5,7), (6,5), (7,4)]

    print(f"The grid number is n = {n}.")

    # --- Step 1: Calculate S_A (Type A Seifert circles) ---
    o_col_to_row = {c: r for c, r in o_pos}
    x_col_to_row = {c: r for c, r in x_pos}
    o_row_to_col = {r: c for c, r in o_pos}

    # pi maps column i to column j where X in col i has same row as O in col j
    pi = {i: o_row_to_col[x_col_to_row[i]] for i in range(1, n + 1)}

    def count_cycles(p, num_elements):
        visited = [False] * (num_elements + 1)
        cycles = 0
        for i in range(1, num_elements + 1):
            if not visited[i]:
                cycles += 1
                j = i
                while not visited[j]:
                    visited[j] = True
                    j = p[j]
        return cycles

    S_A = count_cycles(pi, n)
    print(f"Number of type A Seifert circles (S_A): {S_A}")
    # The result should be 1 for a knot.

    # --- Step 2: Calculate S_B (Type B Seifert circles) ---
    # This involves building a graph where nodes are markers (O_i, X_i)
    # and edges connect markers in the same row or same column.
    x_row_to_col = {r: c for c, r in x_pos}
    
    # Graph has 2n nodes, labeled 1..n for O_i and n+1..2n for X_i
    adj = {}
    for i in range(1, n + 1):
        # Vertical connections: O_i -> X_i
        adj[i] = i + n
        adj[i + n] = i
    
    for r in range(1, n + 1):
        # Horizontal connections: O_i -> X_j where O_i and X_j are in row r
        o_col = o_row_to_col[r]
        x_col = x_row_to_col[r]
        # Connect O_{o_col} and X_{x_col}
        # In our graph, connect node o_col to node x_col+n
        adj[o_col] = x_col + n
        adj[x_col + n] = o_col
    
    # The permutation for this checkerboard graph (nodes 1 to 2n)
    checkerboard_pi = adj
    
    S_B = count_cycles(checkerboard_pi, 2 * n)
    print(f"Number of type B Seifert circles (S_B): {S_B}")

    # --- Step 3: Calculate the number of crossings (c) ---
    row_to_cols_o = {r: c for c, r in o_pos}
    row_to_cols_x = {r: c for c, r in x_pos}
    
    c = 0
    crossing_points = []
    for i in range(1, n + 1):  # column of vertical segment
        for j in range(1, n + 1):  # row of horizontal segment
            o_row_v, x_row_v = o_col_to_row[i], x_col_to_row[i]
            o_col_h, x_col_h = row_to_cols_o[j], row_to_cols_x[j]
            
            is_vertical_span = min(o_row_v, x_row_v) < j < max(o_row_v, x_row_v)
            is_horizontal_span = min(o_col_h, x_col_h) < i < max(o_col_h, x_col_h)
            
            if is_vertical_span and is_horizontal_span:
                c += 1
                crossing_points.append((i,j))

    print(f"Number of crossings (c): {c}")

    # --- Step 4: Calculate the Seifert genus (g) ---
    g_val_times_2 = c - S_A - S_B + 2
    g = g_val_times_2 / 2

    print("\nThe Seifert genus (g) is calculated using the formula: 2g = c - S_A - S_B + 2")
    print(f"2g = {c} - {S_A} - {S_B} + 2 = {g_val_times_2}")
    print(f"So, the genus g = {int(g)}.")

    # --- Step 5: Determine the braid index (b) ---
    print("\nThe braid index (b) is bounded by the genus and the grid number (n):")
    print("2g + 1 <= b <= n")
    
    lower_bound = int(2 * g + 1)
    upper_bound = n
    
    print(f"Substituting the values: 2 * {int(g)} + 1 <= b <= {n}")
    print(f"The inequality simplifies to: {lower_bound} <= b <= {upper_bound}")
    
    if lower_bound == upper_bound:
        print(f"\nTherefore, the braid index of the corresponding knot is {lower_bound}.")
    else:
        print("\nThe bounds are not tight enough to determine the exact braid index with this method.")
    
solve_braid_index()