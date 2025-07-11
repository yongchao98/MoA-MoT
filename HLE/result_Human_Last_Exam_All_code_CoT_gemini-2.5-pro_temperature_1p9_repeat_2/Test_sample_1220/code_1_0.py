def solve_langtons_ant_period():
    """
    Simulates Langton's Ant on a toroidal grid to find its period.

    The period is defined as the number of steps it takes for the ant to return
    to its starting position, facing its starting direction.
    """
    # Grid dimensions
    rows = 4
    cols = 5

    # Initialize a grid of all white squares (0 = white, 1 = black)
    grid = [[0 for _ in range(cols)] for _ in range(rows)]

    # Define the ant's starting state
    start_row, start_col = 0, 0
    # 0: Up, 1: Right, 2: Down, 3: Left
    start_dir = 0  

    # Initialize the ant's current state
    ant_row, ant_col = start_row, start_col
    ant_dir = start_dir

    # Direction vectors for movement [Up, Right, Down, Left]
    # (change in row, change in col)
    d_row = [-1, 0, 1, 0]
    d_col = [0, 1, 0, -1]

    steps = 0
    while True:
        # Check the color of the current square
        if grid[ant_row][ant_col] == 0:  # White square
            # Turn 90° clockwise
            ant_dir = (ant_dir + 1) % 4
            # Flip the color of the square to black
            grid[ant_row][ant_col] = 1
            # Move forward
            ant_row = (ant_row + d_row[ant_dir]) % rows
            ant_col = (ant_col + d_col[ant_dir]) % cols
        else:  # Black square
            # Turn 90° counter-clockwise
            ant_dir = (ant_dir - 1 + 4) % 4
            # Flip the color of the square to white
            grid[ant_row][ant_col] = 0
            # Move forward
            ant_row = (ant_row + d_row[ant_dir]) % rows
            ant_col = (ant_col + d_col[ant_dir]) % cols

        steps += 1

        # Check if the ant has returned to its initial state (position and direction).
        # This condition marks the end of a period for this specific problem.
        if ant_row == start_row and ant_col == start_col and ant_dir == start_dir:
            break

    # The prompt requests to "output each number in the final equation".
    # As there is no equation, we will just print the final result.
    print(steps)

solve_langtons_ant_period()