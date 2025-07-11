def langtons_ant_period():
    """
    Simulates Langton's Ant on a 4x5 toroidal grid to find its period.
    The ant starts at (0, 0) facing up on an all-white grid.
    """
    rows, cols = 4, 5

    # Initialize the grid: 0 for white, 1 for black.
    grid = [[0 for _ in range(cols)] for _ in range(rows)]

    # Initialize the ant's state.
    # Position: [row, col]
    ant_pos = [0, 0]
    # Direction: 0: up, 1: right, 2: down, 3: left
    ant_dir = 0

    # Deltas for position changes corresponding to directions [up, right, down, left]
    dr = [-1, 0, 1, 0]
    dc = [0, 1, 0, -1]

    # History to store states and detect cycles.
    # A state is (ant_row, ant_col, ant_dir, grid_tuple).
    history = {}
    steps = 0

    while True:
        # Create a hashable representation of the grid state.
        grid_tuple = tuple(tuple(row) for row in grid)
        current_state = (ant_pos[0], ant_pos[1], ant_dir, grid_tuple)

        # Check if the state has been seen before.
        if current_state in history:
            period = steps - history[current_state]
            print(period)
            break
        
        # Store the new state and the current step count.
        history[current_state] = steps

        # Get the color of the current square.
        r, c = ant_pos[0], ant_pos[1]
        color = grid[r][c]

        # Apply Langton's Ant rules.
        if color == 0:  # White square
            ant_dir = (ant_dir + 1) % 4  # Turn 90° clockwise
        else:  # Black square
            ant_dir = (ant_dir - 1 + 4) % 4  # Turn 90° counter-clockwise

        # Flip the color of the square.
        grid[r][c] = 1 - color

        # Move the ant forward one unit (with toroidal wrap-around).
        ant_pos[0] = (ant_pos[0] + dr[ant_dir]) % rows
        ant_pos[1] = (ant_pos[1] + dc[ant_dir]) % cols
        
        steps += 1

langtons_ant_period()