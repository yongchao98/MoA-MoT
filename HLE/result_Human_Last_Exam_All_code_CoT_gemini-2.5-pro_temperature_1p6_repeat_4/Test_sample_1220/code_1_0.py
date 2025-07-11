def solve_langtons_ant_period():
    """
    Simulates Langton's Ant on a toroidal grid to find its period.

    The period is the number of steps required for the system (ant's position,
    ant's direction, and grid colors) to return to its initial state.
    """
    # 1. Initialization
    rows = 4
    cols = 5

    # Ant's state: position (x, y) and direction.
    # Direction: 0:Up, 1:Right, 2:Down, 3:Left
    ant_x, ant_y = 0, 0
    direction = 0  # Starts facing Up

    # Store the initial state for period detection
    initial_ant_x, initial_ant_y, initial_direction = ant_x, ant_y, direction

    # Grid: 2D list where 0 is white and 1 is black.
    grid = [[0 for _ in range(cols)] for _ in range(rows)]
    # A deep copy of the initial grid to compare against for period detection.
    initial_grid = [row[:] for row in grid]

    # Movement vectors (dx, dy) for directions [Up, Right, Down, Left]
    # Grid coordinates are (row, col), which corresponds to (y, x).
    # Up: (y-1), Right: (x+1), Down: (y+1), Left: (x-1)
    move_x = [0, 1, 0, -1]
    move_y = [-1, 0, 1, 0]

    step_count = 0

    # 2. Simulation Loop
    while True:
        # Rules are applied based on the color of the current square
        is_white_square = (grid[ant_y][ant_x] == 0)

        if is_white_square:
            # At a white square: turn 90° clockwise, flip color, move forward
            direction = (direction + 1) % 4
            grid[ant_y][ant_x] = 1
        else:  # Black square
            # At a black square: turn 90° counter-clockwise, flip color, move forward
            direction = (direction - 1 + 4) % 4
            grid[ant_y][ant_x] = 0

        # Move the ant one unit forward in its new direction
        ant_x = (ant_x + move_x[direction]) % cols
        ant_y = (ant_y + move_y[direction]) % rows

        step_count += 1

        # 3. Check for return to the complete initial state
        if (ant_x == initial_ant_x and
            ant_y == initial_ant_y and
            direction == initial_direction and
            grid == initial_grid):
            break

    # 4. Output the Result
    print(f"The period for Langton's Ant on a {rows} by {cols} torus is calculated.")
    print(f"Final calculation: Period = {step_count}")

solve_langtons_ant_period()
<<<604>>>