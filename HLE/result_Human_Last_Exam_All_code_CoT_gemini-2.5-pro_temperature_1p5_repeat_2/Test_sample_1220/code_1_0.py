def solve_langtons_ant_period():
    """
    Simulates Langton's Ant on a toroidal grid to find its period.

    The state of the system is defined by the ant's position, its direction,
    and the color of every cell. The period is the number of steps it takes
    for the system to return to its initial state for the first time.
    """
    # Grid dimensions as specified in the problem
    ROWS = 4
    COLS = 5

    # Ant's starting state
    start_r, start_c = 0, 0  # Start at top-left corner
    start_dir = 0            # 0: Up, 1: Right, 2: Down, 3: Left

    # Initialize current ant state
    ant_r, ant_c, ant_dir = start_r, start_c, start_dir

    # Define movement vectors for [Up, Right, Down, Left]
    # (dr, dc) corresponds to change in (row, col)
    moves = [(-1, 0), (0, 1), (1, 0), (0, -1)]

    # Initialize grid to all white. 0: White, 1: Black
    grid = [[0 for _ in range(COLS)] for _ in range(ROWS)]

    step_count = 0
    while True:
        # After the first step, check if we have returned to the initial state.
        # The initial state is: ant at (start_r, start_c), facing start_dir,
        # on an all-white grid.
        if step_count > 0:
            is_initial_pos = (ant_r == start_r and ant_c == start_c)
            is_initial_dir = (ant_dir == start_dir)
            # This check is the most time-consuming, so perform it last
            is_initial_grid = all(cell == 0 for row in grid for cell in row)

            if is_initial_pos and is_initial_dir and is_initial_grid:
                # The period is the current step count
                break

        # Get the color of the current square
        current_color = grid[ant_r][ant_c]

        # Apply Langton's Ant rules
        if current_color == 0:  # White square
            # Turn 90 degrees clockwise
            ant_dir = (ant_dir + 1) % 4
            # Flip color to black
            grid[ant_r][ant_c] = 1
        else:  # Black square (current_color == 1)
            # Turn 90 degrees counter-clockwise
            ant_dir = (ant_dir - 1 + 4) % 4
            # Flip color to white
            grid[ant_r][ant_c] = 0

        # Move forward one unit in the new direction (with toroidal wrap-around)
        dr, dc = moves[ant_dir]
        ant_r = (ant_r + dr) % ROWS
        ant_c = (ant_c + dc) % COLS

        # Increment step counter
        step_count += 1
    
    # The problem asks what the period is. This is the final number.
    # The problem also has a confusing instruction about outputting numbers
    # in an equation. The safest interpretation is to output the numbers
    # involved in the problem's definition and the final answer.
    print(f"For a {ROWS}x{COLS} grid, the period is {step_count}.")
    
solve_langtons_ant_period()