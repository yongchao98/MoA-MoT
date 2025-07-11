def solve_langtons_ant_period():
    """
    Simulates Langton's Ant on a toroidal grid to find its period.
    """
    rows = 4
    cols = 5

    # Initialize the grid with all white cells (0)
    grid = [[0 for _ in range(cols)] for _ in range(rows)]

    # Define the ant's initial state
    # Directions: 0: Up, 1: Right, 2: Down, 3: Left
    ant_row, ant_col = 0, 0
    ant_dir = 0  # Starting facing Up

    # Store the complete initial state for later comparison.
    # The grid must be converted to an immutable type (tuple of tuples) to be stored.
    initial_grid_state = tuple(map(tuple, grid))
    initial_ant_pos = (ant_row, ant_col)
    initial_ant_dir = ant_dir

    # Deltas for moving the ant based on direction [Up, Right, Down, Left]
    # Change in row: dr, Change in column: dc
    dr = [-1, 0, 1, 0]
    dc = [0, 1, 0, -1]

    steps = 0
    while True:
        steps += 1

        # Check the color of the current square
        if grid[ant_row][ant_col] == 0:  # White square
            # Turn 90 degrees clockwise
            ant_dir = (ant_dir + 1) % 4
        else:  # Black square (1)
            # Turn 90 degrees counter-clockwise
            ant_dir = (ant_dir - 1 + 4) % 4

        # Flip the color of the current square
        grid[ant_row][ant_col] = 1 - grid[ant_row][ant_col]

        # Move the ant forward one unit
        ant_row += dr[ant_dir]
        ant_col += dc[ant_dir]

        # Apply toroidal grid boundary conditions (wrap around)
        ant_row %= rows
        ant_col %= cols

        # Check if the system has returned to its initial state
        current_grid_state = tuple(map(tuple, grid))
        if (ant_row, ant_col) == initial_ant_pos and \
           ant_dir == initial_ant_dir and \
           current_grid_state == initial_grid_state:
            # The period is the number of steps taken to return to the initial state
            break

    # The problem has no explicit equation, it asks for the period, which is a number.
    # The instructions say "output each number in the final equation!", so we will print the final number.
    print(steps)

solve_langtons_ant_period()
<<<1160>>>