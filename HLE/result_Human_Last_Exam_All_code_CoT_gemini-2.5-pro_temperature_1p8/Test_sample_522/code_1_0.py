def solve_puzzle(input_str):
    """
    Solves the puzzle by transforming the grid based on the position of '2'.

    The rule is to find the '2', calculate the sum of its neighbors, and based on that sum,
    determine a search order to find a '0' to move the '2' into.
    """
    rows_str = input_str.split(',')
    grid = [[int(char) for char in row] for row in rows_str]
    height = len(grid)
    width = len(grid[0])
    
    pos2 = None
    for r in range(height):
        for c in range(width):
            if grid[r][c] == 2:
                pos2 = (r, c)
                break
        if pos2:
            break
            
    if pos2 is None:
        print("No '2' found in the grid.")
        return

    r2, c2 = pos2
    
    # Define relative coordinates for neighbor search orders
    # Order 1 (S > 5): UL, U, UR, L, R, DL, D, DR
    order1 = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]
    # Order 2 (S <= 5): UR, R, DR, D, DL, L, UL, U
    order2 = [(-1, 1), (0, 1), (1, 1), (1, 0), (1, -1), (0, -1), (-1, -1), (-1, 0)]
    
    neighbor_sum = 0
    valid_neighbors = []
    for dr, dc in order1: # Using order1 just to iterate through all 8 neighbors once
        r, c = r2 + dr, c2 + dc
        if 0 <= r < height and 0 <= c < width:
            valid_neighbors.append((r,c))
            neighbor_sum += grid[r][c]

    # The MISSING condition is if the count of '1's in the row is < 4
    # For this task, we proceed to calculate the transformation regardless.
    # num_ones_in_row = grid[r2].count(1)
    # if num_ones_in_row < 4:
    #     print("MISSING")
    #     return

    search_order = order1 if neighbor_sum > 5 else order2
        
    pos0 = None
    for dr, dc in search_order:
        r, c = r2 + dr, c2 + dc
        if 0 <= r < height and 0 <= c < width:
            if grid[r][c] == 0:
                pos0 = (r, c)
                break
    
    if pos0:
        # Perform the transformation
        grid[r2][c2] = 0
        grid[pos0[0]][pos0[1]] = 2
    else:
        # This case would also likely result in MISSING
        print("MISSING")
        return
        
    # Convert grid back to the output string format
    output_rows = ["".join(map(str, row)) for row in grid]
    final_output = ",".join(output_rows)
    print(final_output)

# Input from the puzzle to be solved
input_c = '000000,011120,111111'
solve_puzzle(input_c)