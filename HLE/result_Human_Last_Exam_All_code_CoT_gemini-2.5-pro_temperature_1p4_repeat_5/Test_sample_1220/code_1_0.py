import sys

def find_langtons_ant_period():
    """
    This script simulates Langton's Ant on a toroidal grid to find its period.

    The simulation proceeds as follows:
    1. A grid of a specified size is initialized with all white cells.
    2. The ant is placed on the grid at a starting position, facing a starting direction.
    3. The initial state of the system (ant's position, ant's direction, and grid colors) is stored.
    4. The simulation runs step by step, applying Langton's Ant rules.
    5. After each step, the current state is compared to the initial state.
    6. The simulation stops when the system returns to its initial state. The number of steps
       taken is the period of the ant's path.
    """
    
    # 1. Define grid dimensions and ant properties
    rows, cols = 4, 5
    
    # Grid initialization: 0 for white, 1 for black
    grid = [[0 for _ in range(cols)] for _ in range(rows)]
    
    # Ant's initial state
    ant_r, ant_c = 0, 0  # Starting at top-left corner
    
    # Directions are represented numerically: 0:Up, 1:Right, 2:Down, 3:Left
    ant_dir = 0  # Starts facing Up as per the problem description
    
    # Deltas for movement based on direction [Up, Right, Down, Left]
    dr = [-1, 0, 1, 0]
    dc = [0, 1, 0, -1]
    
    # 2. Store the initial state to detect when a full cycle is complete.
    # The state must include the ant's position, direction, and the entire grid's configuration.
    # We convert the grid (a list of lists) to a tuple of tuples to make it hashable and comparable.
    def get_current_state():
        return ((ant_r, ant_c, ant_dir), tuple(map(tuple, grid)))

    initial_state = get_current_state()
    steps = 0
    
    # 3. Start the simulation loop
    while True:
        # Get the color of the current square where the ant is located
        current_color = grid[ant_r][ant_c]
        
        # Apply Langton's Ant rules based on the color
        if current_color == 0:  # If the square is white
            # Turn 90 degrees clockwise
            ant_dir = (ant_dir + 1) % 4
        else:  # If the square is black
            # Turn 90 degrees counter-clockwise
            ant_dir = (ant_dir - 1 + 4) % 4
            
        # Flip the color of the current square
        grid[ant_r][ant_c] = 1 - current_color
        
        # Move the ant forward one unit in its new direction
        # The modulo operator (%) ensures the grid is toroidal (wraps around)
        ant_r = (ant_r + dr[ant_dir]) % rows
        ant_c = (ant_c + dc[ant_dir]) % cols
        
        # Increment the step counter
        steps += 1
        
        # Check if the system has returned to its initial state
        if get_current_state() == initial_state:
            break # Period found, exit the loop
            
        # A safety break to prevent an infinite loop in case of an error.
        # This value is set high enough not to interfere with the expected result.
        if steps > 500000:
            print("Error: Simulation exceeded 500,000 steps.", file=sys.stderr)
            return

    # 4. Print the final result
    print("The period of Langton's Ant on a 4x5 toroidal grid is:")
    print(steps)

if __name__ == '__main__':
    find_langtons_ant_period()
<<<264>>>