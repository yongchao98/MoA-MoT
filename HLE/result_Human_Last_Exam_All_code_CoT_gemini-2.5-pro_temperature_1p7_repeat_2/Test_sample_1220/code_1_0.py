import sys

def solve_langtons_ant_period():
    """
    Simulates Langton's Ant on a 4x5 toroidal grid to find the period of its path.

    The simulation tracks the complete state of the system (ant's position, ant's
    direction, and the grid's color pattern) at each step. It uses a dictionary
    to store the history of these states. When a state repeats, a cycle is
    detected, and the period is the length of that cycle.
    """
    rows, cols = 4, 5

    # --- State Initialization ---
    # The grid is represented by a list of lists. 0 is white, 1 is black.
    grid = [[0 for _ in range(cols)] for _ in range(rows)]
    
    # The ant starts at the top-left corner.
    ant_row, ant_col = 0, 0
    
    # The ant's direction is represented by an integer.
    # 0: up, 1: right, 2: down, 3: left
    ant_dir = 0  # Starts facing "up"

    # --- Cycle Detection ---
    # `history` stores visited states. Key: state tuple, Value: step number.
    history = {}
    step = 0

    # --- Simulation Constants ---
    # Directional vectors for movement (change in row, change in col).
    dr = [-1, 0, 1, 0]  # Corresponds to up, right, down, left
    dc = [0, 1, 0, -1]  # Corresponds to up, right, down, left

    while True:
        # Create a hashable representation of the current state.
        # The grid list must be converted to a tuple of tuples to be used as a key.
        grid_tuple = tuple(map(tuple, grid))
        current_state = (ant_row, ant_col, ant_dir, grid_tuple)

        # Check if this exact state has been seen before.
        if current_state in history:
            start_of_cycle_step = history[current_state]
            period = step - start_of_cycle_step

            # --- Output Results ---
            print("A repeating state has been found, meaning the ant's path is periodic.")
            print(f"The state at step {start_of_cycle_step} is identical to the state at step {step}.")
            print("The period of the ant's path is the length of this cycle, calculated as:")
            print(f"{step} - {start_of_cycle_step} = {period}")
            return period

        # If the state is new, record it in the history.
        history[current_state] = step

        # --- Apply Langton's Ant Rules ---
        current_color = grid[ant_row][ant_col]

        # 1. Flip the color of the current square.
        grid[ant_row][ant_col] = 1 - current_color

        # 2. Turn the ant based on the square's original color.
        if current_color == 0:  # White square
            ant_dir = (ant_dir + 1) % 4  # Turn 90 degrees clockwise
        else:  # Black square
            ant_dir = (ant_dir - 1 + 4) % 4  # Turn 90 degrees counter-clockwise

        # 3. Move the ant forward one unit (with toroidal wrap-around).
        ant_row = (ant_row + dr[ant_dir]) % rows
        ant_col = (ant_col + dc[ant_dir]) % cols
        
        step += 1

if __name__ == '__main__':
    solve_langtons_ant_period()