def solve_puzzle(input_str):
    """
    Solves the puzzle by moving the '2' piece diagonally based on a set of rules.
    """
    # Parse the input string into a 2D list (grid)
    grid = [list(row) for row in input_str.split(',')]
    rows = len(grid)
    if rows == 0:
        print(input_str)
        return
    cols = len(grid[0])

    # Find the initial position of '2'
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

    # If no '2' is found, no changes are made
    if not found_2:
        print(input_str)
        return

    # Define the move preference order: up-right, up-left, down-right, down-left
    moves = [
        (-1, 1),   # Up-Right
        (-1, -1),  # Up-Left
        (1, 1),    # Down-Right
        (1, -1)    # Down-Left
    ]

    # Find and execute the first valid move
    move_made = False
    for dr, dc in moves:
        nr, nc = r2 + dr, c2 + dc

        # Check if the move is within bounds
        if 0 <= nr < rows and 0 <= nc < cols:
            # Check if the destination is a '0'
            if grid[nr][nc] == '0':
                # Perform the move
                grid[r2][c2] = '0'
                grid[nr][nc] = '2'
                move_made = True
                break  # Exit after the first successful move

    # Convert the final grid back to a string and print it
    output_rows = ["".join(row) for row in grid]
    result = ",".join(output_rows)
    print(result)

# The input for which we need to find the MISSING value
input_c = '000000,011120,111111'
solve_puzzle(input_c)