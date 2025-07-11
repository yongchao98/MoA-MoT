def solve_puzzle():
    """
    Solves the puzzle by finding and applying the transformation rule.
    The rule involves moving a '2' to an adjacent '0' based on a conditional
    directional preference determined by the '2's neighbors.
    """
    input_str = '000000,011120,111111'

    # Parse the input string into a 2D list of characters (grid)
    grid = [list(row) for row in input_str.split(',')]
    rows = len(grid)
    cols = len(grid[0])

    # Find the starting coordinates of the '2'
    start_r, start_c = -1, -1
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == '2':
                start_r, start_c = r, c
                break
        if start_r != -1:
            break
    
    # This case is not expected based on the problem description
    if start_r == -1:
        print(f"Error: No '2' found in the input grid.")
        return

    # Find all '0'-valued neighbors and check the "East flank"
    neighbor_offsets = [
        (-1, -1), (-1, 0), (-1, 1),  # Top row (NW, N, NE)
        (0, -1),           (0, 1),   # Middle row (W, E)
        (1, -1), (1, 0), (1, 1)    # Bottom row (SW, S, SE)
    ]
    
    zero_neighbors = []
    east_flank_has_zero = False

    for dr, dc in neighbor_offsets:
        nr, nc = start_r + dr, start_c + dc
        
        # Check if the neighbor is within the grid boundaries
        if 0 <= nr < rows and 0 <= nc < cols:
            if grid[nr][nc] == '0':
                zero_neighbors.append((nr, nc))
                # Check if this neighbor is on the "East flank" (column offset is +1)
                if dc == 1:
                    east_flank_has_zero = True
    
    if not zero_neighbors:
        # This case is not expected for the given puzzle, as a move seems implied.
        # If it occurred, the grid would not change.
        final_result = ",".join(["".join(row) for row in grid])
        print(final_result)
        return

    # Determine the destination based on the conditional rule
    dest_r, dest_c = -1, -1
    # Find the minimum row index among all '0'-valued neighbors
    min_r = min(r for r, c in zero_neighbors)
    # Filter candidates to only those in the top-most row
    top_row_candidates = [c for r, c in zero_neighbors if r == min_r]
    
    if east_flank_has_zero:
        # "Rightward" rule: minimum row, maximum column
        dest_c = max(top_row_candidates)
    else:
        # "Leftward" rule: minimum row, minimum column
        dest_c = min(top_row_candidates)
    
    dest_r = min_r

    # Create the new grid by applying the transformation
    grid[start_r][start_c] = '0'
    grid[dest_r][dest_c] = '2'

    # Join the grid back into the final string format and print
    output_rows = ["".join(row) for row in grid]
    final_result = ",".join(output_rows)
    print(final_result)

solve_puzzle()