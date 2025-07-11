def find_langtons_ant_period():
    """
    Simulates Langton's Ant on a 4x5 toroidal grid to find the period of its path.

    The state of the system is defined by the ant's position, the ant's direction,
    and the color of every cell on the grid. The simulation runs until a state
    repeats, which indicates a cycle has been found. The period is the length
    of this cycle.
    """
    rows, cols = 4, 5

    # Ant's state: position [row, col] and direction.
    # Directions are encoded as: 0: Up, 1: Right, 2: Down, 3: Left.
    ant_pos = [0, 0]  # Start at top-left corner
    ant_dir = 0       # Start facing Up

    # The grid is initialized to all white squares (represented by 0).
    # Black squares are represented by 1.
    grid = [[0 for _ in range(cols)] for _ in range(rows)]

    # 'history' stores encountered states to detect a cycle.
    # Key: A tuple representing the state (ant_row, ant_col, ant_dir, grid_tuple)
    # Value: The step number when the state was first seen.
    history = {}
    steps = 0

    while True:
        # Create a hashable (and thus storable in a dictionary) representation of the grid.
        grid_tuple = tuple(tuple(row) for row in grid)

        # The complete state of the system.
        current_state = (ant_pos[0], ant_pos[1], ant_dir, grid_tuple)

        # If the state has been seen before, a cycle is found.
        if current_state in history:
            first_occurrence_step = history[current_state]
            period = steps - first_occurrence_step
            # Print the final equation for the period calculation.
            print(f"The ant's path enters a cycle of length {period}.")
            print(f"Cycle detected at step {steps}, which previously occurred at step {first_occurrence_step}.")
            print(f"Period calculation: {steps} - {first_occurrence_step} = {period}")
            return period

        # Store the current state and step number.
        history[current_state] = steps

        # --- Apply Langton's Ant rules ---
        
        # 1. Check color and turn accordingly.
        current_color = grid[ant_pos[0]][ant_pos[1]]
        if current_color == 0:  # White square
            ant_dir = (ant_dir + 1) % 4  # Turn 90 degrees clockwise
        else:  # Black square
            ant_dir = (ant_dir - 1 + 4) % 4  # Turn 90 degrees counter-clockwise

        # 2. Flip the color of the current square.
        grid[ant_pos[0]][ant_pos[1]] = 1 - current_color

        # 3. Move the ant one step forward.
        if ant_dir == 0:  # Up
            ant_pos[0] = (ant_pos[0] - 1 + rows) % rows
        elif ant_dir == 1:  # Right
            ant_pos[1] = (ant_pos[1] + 1) % cols
        elif ant_dir == 2:  # Down
            ant_pos[0] = (ant_pos[0] + 1) % rows
        elif ant_dir == 3:  # Left
            ant_pos[1] = (ant_pos[1] - 1 + cols) % cols
            
        steps += 1

# Execute the simulation and find the period.
period = find_langtons_ant_period()
<<<188>>>