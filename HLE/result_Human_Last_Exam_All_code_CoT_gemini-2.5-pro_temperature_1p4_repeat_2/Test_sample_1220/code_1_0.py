import sys

def solve_langtons_ant_period():
    """
    Simulates Langton's Ant on a toroidal grid to find its period.

    The state of the system is defined by the ant's position, its direction,
    and the color of every cell on the grid. The simulation runs until a
    state repeats, and the number of steps taken is the period.
    """
    # Grid dimensions
    ROWS = 4
    COLS = 5

    # Ant's initial state
    ant_r, ant_c = 0, 0  # Start at top-left corner
    # Directions: 0:Up, 1:Right, 2:Down, 3:Left
    ant_dir = 0  # Starts facing Up

    # Grid initialization: 0 for white, 1 for black. Starts all white.
    grid = [[0 for _ in range(COLS)] for _ in range(ROWS)]

    # Direction vectors for movement [dr, dc]
    # Corresponds to directions 0, 1, 2, 3
    moves = [(-1, 0), (0, 1), (1, 0), (0, -1)]

    # History of states to detect a cycle.
    # Key: state tuple (ant_r, ant_c, ant_dir, grid_tuple)
    # Value: step number
    history = {}
    steps = 0

    while True:
        # Create a hashable representation of the grid for the state key
        grid_tuple = tuple(tuple(row) for row in grid)
        current_state = (ant_r, ant_c, ant_dir, grid_tuple)

        # Check if we have returned to a previously seen state
        if current_state in history:
            period = steps
            break

        # Record the new state and the current step count
        history[current_state] = steps

        # Get the color of the ant's current square
        color = grid[ant_r][ant_c]

        # Apply Langton's Ant rules
        if color == 0:  # White square
            # Turn 90 degrees clockwise
            ant_dir = (ant_dir + 1) % 4
            # Flip the color of the square to black
            grid[ant_r][ant_c] = 1
        else:  # Black square (color == 1)
            # Turn 90 degrees counter-clockwise
            ant_dir = (ant_dir - 1 + 4) % 4
            # Flip the color of the square to white
            grid[ant_r][ant_c] = 0

        # Move forward one unit based on the new direction
        dr, dc = moves[ant_dir]
        # Apply toroidal wrapping
        ant_r = (ant_r + dr) % ROWS
        ant_c = (ant_c + dc) % COLS

        # Increment the step counter
        steps += 1
        
        # A safety break for very long running simulations, though not expected for this size.
        if steps > 20000:
            print("Simulation exceeded maximum steps.", file=sys.stderr)
            return

    print(f"The period of Langton's Ant on a {ROWS} by {COLS} toroidal grid is: {period}")

solve_langtons_ant_period()
<<<1240>>>