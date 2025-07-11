def solve_puzzle():
    """
    Solves the puzzle by applying a transformation rule to a grid.

    The rule involves an agent '2' moving to an adjacent '0' cell based on a
    priority list of directions. The output is 'MISSING' if the move is
    ambiguous (the top two priority moves are both possible) or impossible.
    """
    # The input for which we need to find the output.
    input_str = '000000,011120,111111'

    # Step 1: Parse the input string into a 2D grid of integers.
    try:
        grid = [[int(char) for char in row] for row in input_str.split(',')]
        height = len(grid)
        if height == 0:
            print("MISSING")
            return
        width = len(grid[0])
    except (ValueError, IndexError):
        # Handle malformed input string
        print("MISSING")
        return

    # Step 2: Find the starting location of the agent '2'.
    start_pos = None
    for r in range(height):
        for c in range(width):
            if grid[r][c] == 2:
                start_pos = (r, c)
                break
        if start_pos:
            break

    # If no '2' is found, the transformation is not applicable.
    if not start_pos:
        print("MISSING")
        return

    r_start, c_start = start_pos

    # Step 3: Define the priority order for movement directions.
    # Order: Top-Right, Right, Bottom-Right, Bottom, Bottom-Left, Left, Top-Left, Top
    directions = [
        (-1, 1), (0, 1), (1, 1), (1, 0),
        (1, -1), (0, -1), (-1, -1), (-1, 0)
    ]

    # Step 4: Find all valid moves.
    # A move is valid if the destination is within bounds and is a '0'.
    valid_moves = []
    for dr, dc in directions:
        nr, nc = r_start + dr, c_start + dc
        if 0 <= nr < height and 0 <= nc < width and grid[nr][nc] == 0:
            valid_moves.append((dr, dc))

    # Step 5: Check for conditions that result in a "MISSING" output.
    # Condition 1: No possible moves (agent is trapped).
    if not valid_moves:
        print("MISSING")
        return

    # Condition 2: Ambiguity. The two highest-priority moves are both possible.
    # Convert valid_moves to a set for efficient lookup.
    valid_moves_set = set(valid_moves)
    if directions[0] in valid_moves_set and directions[1] in valid_moves_set:
        print("MISSING")
        return

    # Step 6: If no 'MISSING' condition is met, perform the highest-priority move.
    # The move to make is the first one from the priority list that is also in valid_moves.
    best_move = None
    for move in directions:
        if move in valid_moves_set:
            best_move = move
            break
            
    dr, dc = best_move
    nr, nc = r_start + dr, c_start + dc
    
    # Create a new grid for the output by applying the move.
    output_grid = [row[:] for row in grid]
    output_grid[r_start][c_start] = 0
    output_grid[nr][nc] = 2

    # Step 7: Format the resulting grid back into a string and print it.
    output_rows = ["".join(map(str, row)) for row in output_grid]
    print(",".join(output_rows))

# Execute the function to find the solution for Input_C.
solve_puzzle()