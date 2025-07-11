import collections

def solve_langtons_ant_period():
    """
    Simulates Langton's Ant on a toroidal grid to find its period.

    The state of the system is defined by the ant's position, its direction,
    and the configuration of the entire grid. The simulation runs step by step,
    storing each unique state encountered. The period is the number of steps
    in the cycle, which is found when a state repeats.
    """
    rows, cols = 4, 5
    grid = [[0 for _ in range(cols)] for _ in range(rows)]

    # Ant's initial state: position (0, 0), direction Up
    ant_row, ant_col = 0, 0
    # Directions: 0: Up, 1: Right, 2: Down, 3: Left
    ant_dir = 0
    
    # Moves corresponding to directions: (d_row, d_col)
    moves = {
        0: (-1, 0),  # Up
        1: (0, 1),   # Right
        2: (1, 0),   # Down
        3: (0, -1)   # Left
    }

    # History to store states and detect a cycle
    # Key: state tuple, Value: step number
    history = {}
    steps = 0

    while True:
        # Create a hashable representation of the grid for the state tuple
        grid_tuple = tuple(tuple(row) for row in grid)
        
        # Define the current state of the entire system
        current_state = (ant_row, ant_col, ant_dir, grid_tuple)

        # Check if this state has been seen before
        if current_state in history:
            # A cycle is found. The period is the number of steps taken since
            # the first time we saw this state.
            period = steps - history[current_state]
            print(period)
            return

        # Store the current state and the step count
        history[current_state] = steps
        
        # Get the color of the current square
        color = grid[ant_row][ant_col]

        if color == 0:  # White square
            # Turn 90 degrees clockwise
            ant_dir = (ant_dir + 1) % 4
            # Flip the color of the square to black
            grid[ant_row][ant_col] = 1
        else:  # Black square (color == 1)
            # Turn 90 degrees counter-clockwise
            ant_dir = (ant_dir - 1 + 4) % 4
            # Flip the color of the square to white
            grid[ant_row][ant_col] = 0

        # Move the ant forward one unit
        move = moves[ant_dir]
        ant_row = (ant_row + move[0]) % rows
        ant_col = (ant_col + move[1]) % cols
            
        # Increment the step counter
        steps += 1

solve_langtons_ant_period()
<<<292>>>