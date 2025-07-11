def solve_puzzle():
    """
    Solves the puzzle by moving a '2' to a diagonal '0' based on a specific search order.
    The rule is: The '2' moves to the first available '0' found by checking its
    diagonal neighbors in the order: Bottom-Right, Bottom-Left, Top-Right, Top-Left.
    The original position of the '2' becomes '0'.
    """
    input_str = '000000,011120,111111'

    # Parse the input string into a 2D list (grid) of characters.
    grid_rows = input_str.split(',')
    grid = [list(row) for row in grid_rows]
    
    rows = len(grid)
    if rows == 0:
        print(input_str)
        return
    cols = len(grid[0])

    # Find the coordinates (r, c) of the '2'.
    r2, c2 = -1, -1
    found_2 = False
    for r in range(rows):
        for c in range(cols):
            if grid[r][c] == '2':
                r2, c2 = r, c
                found_2 = True
                break
        if found_2:
            break

    # If no '2' is found, the grid remains unchanged.
    if not found_2:
        print(input_str)
        return

    # Define the diagonal move search order: (dr, dc)
    # 1. Bottom-Right, 2. Bottom-Left, 3. Top-Right, 4. Top-Left
    moves = [(1, 1), (1, -1), (-1, 1), (-1, -1)]

    # Find the first valid move and apply it.
    for dr, dc in moves:
        new_r, new_c = r2 + dr, c2 + dc

        # Check if the new coordinates are within the grid boundaries.
        if 0 <= new_r < rows and 0 <= new_c < cols:
            # Check if the destination cell contains a '0'.
            if grid[new_r][new_c] == '0':
                # Perform the move.
                grid[r2][c2] = '0'
                grid[new_r][new_c] = '2'
                # Stop after the first successful move.
                break
    
    # Convert the modified grid back to the comma-separated string format.
    output_rows = ["".join(row) for row in grid]
    result = ",".join(output_rows)
    
    print(result)

solve_puzzle()