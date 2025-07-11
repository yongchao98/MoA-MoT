import sys

def solve_langtons_ant_period():
    """
    Simulates Langton's Ant on a toroidal grid to find its period.

    The period is the number of steps it takes for the entire system (ant position,
    ant direction, and grid configuration) to return to its initial state.
    """
    # 1. Initialize the system state
    rows, cols = 4, 5
    # The grid is initialized to all white (0). Black is represented by 1.
    grid = [[0] * cols for _ in range(rows)]

    # The ant starts at the top-left corner, (0, 0). The choice is arbitrary
    # as the problem statement says.
    ant_row, ant_col = 0, 0

    # Directions are represented numerically: 0: Up, 1: Right, 2: Down, 3: Left.
    # The ant starts facing 'up'.
    ant_dir = 0
    # The change in (row, col) corresponding to each direction.
    moves = [(-1, 0), (0, 1), (1, 0), (0, -1)]  # Up, Right, Down, Left

    # We must store the initial state to know when the cycle completes.
    # The grid configuration is converted to a hashable tuple to be stored.
    initial_ant_pos = (ant_row, ant_col)
    initial_ant_dir = ant_dir
    initial_grid_tuple = tuple(tuple(row) for row in grid)

    steps = 0

    # 2. Run the simulation loop
    while True:
        # At the start of each loop, check if we have returned to the initial state.
        # This check is skipped at step 0. The number of steps taken is the period.
        if steps > 0:
            # We must check all three components of the state.
            if (ant_row, ant_col) == initial_ant_pos and \
               ant_dir == initial_ant_dir and \
               tuple(tuple(row) for row in grid) == initial_grid_tuple:
                # The system has returned to its initial state. The period is found.
                break

        # Check for excessively long runtimes as a safeguard.
        if steps > 2000:
             print("Error: Simulation is taking longer than expected.", file=sys.stderr)
             return -1

        # 3. Apply Langton's Ant rules
        current_color = grid[ant_row][ant_col]

        if current_color == 0:  # If the square is white
            # Turn 90 degrees clockwise.
            ant_dir = (ant_dir + 1) % 4
            # Flip the color of the square to black.
            grid[ant_row][ant_col] = 1
        else:  # If the square is black (1)
            # Turn 90 degrees counter-clockwise.
            ant_dir = (ant_dir - 1 + 4) % 4
            # Flip the color of the square to white.
            grid[ant_row][ant_col] = 0

        # Move the ant forward one unit in its new direction.
        dr, dc = moves[ant_dir]
        
        # The modulo operator (%) handles the toroidal nature of the grid,
        # wrapping the ant around to the opposite edge if it moves off.
        ant_row = (ant_row + dr) % rows
        ant_col = (ant_col + dc) % cols

        # Increment the step counter.
        steps += 1
        
    return steps

# Calculate the period by running the simulation.
rows = 4
cols = 5
period = solve_langtons_ant_period()

if period != -1:
    print(f"The period of the ant on a torus with {rows} rows and {cols} columns is: {period}")

<<<608>>>