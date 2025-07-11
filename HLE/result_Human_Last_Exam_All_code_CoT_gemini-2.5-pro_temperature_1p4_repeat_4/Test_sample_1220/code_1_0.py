def find_langtons_ant_period():
    """
    This function simulates Langton's Ant on a 4x5 toroidal grid to find its period.

    The state of the simulation is defined by the ant's position, its direction,
    and the configuration of the grid. The simulation runs step by step, recording
    each unique state in a history dictionary. When a state repeats, a cycle is
    detected. The period of this cycle is the difference between the current step
    and the step number of the state's first occurrence.
    """
    rows, cols = 4, 5

    # 0: white, 1: black
    grid = [[0 for _ in range(cols)] for _ in range(rows)]

    # Ant's initial state: position (0,0), direction Up
    ant_row, ant_col = 0, 0
    # 0: Up, 1: Right, 2: Down, 3: Left
    ant_dir = 0

    #           Up      Right    Down    Left
    moves = [(-1, 0), (0, 1), (1, 0), (0, -1)]

    # History stores {state: step_number}
    history = {}
    steps = 0

    while True:
        # Create a hashable tuple representation of the grid for use as a dictionary key
        grid_tuple = tuple(tuple(row) for row in grid)
        current_state = (ant_row, ant_col, ant_dir, grid_tuple)

        # Check if we have been in this state before
        if current_state in history:
            first_occurrence_step = history[current_state]
            period = steps - first_occurrence_step
            print("Cycle Detected.")
            print(f"The state at step {steps} is a repeat of the state from step {first_occurrence_step}.")
            print(f"The period is the current step minus the step of the first occurrence: {steps} - {first_occurrence_step} = {period}")
            print(period)
            break

        # Record the new state
        history[current_state] = steps

        # Get the color of the current square
        color = grid[ant_row][ant_col]

        # Apply Langton's Ant rules
        if color == 0:  # White square
            ant_dir = (ant_dir + 1) % 4  # Turn 90° clockwise
            grid[ant_row][ant_col] = 1   # Flip color to black
        else:  # Black square
            ant_dir = (ant_dir - 1 + 4) % 4  # Turn 90° counter-clockwise
            grid[ant_row][ant_col] = 0       # Flip color to white

        # Move forward one unit based on the new direction
        dr, dc = moves[ant_dir]
        ant_row += dr
        ant_col += dc

        # Apply toroidal grid logic (wrap around)
        ant_row %= rows
        ant_col %= cols

        # Increment step counter
        steps += 1

if __name__ == '__main__':
    find_langtons_ant_period()
<<<264>>>