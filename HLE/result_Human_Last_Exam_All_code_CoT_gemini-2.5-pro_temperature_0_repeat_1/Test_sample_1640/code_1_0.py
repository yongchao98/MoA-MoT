import copy

def solve_tetris_puzzle():
    """
    This function simulates a Tetris-like game to find which command sequence
    results in a specific target grid configuration.
    """
    
    # --- Game Configuration ---
    GRID_WIDTH = 10
    GRID_HEIGHT = 10
    START_COL = 4
    EMPTY_CELL_VAL = 0
    FILLED_CELL_VAL = 1

    # Tetromino shapes are defined as matrices of 0s and 1s.
    TETROMINOES = {
        'straight': [[1, 1, 1, 1]],
        'square': [[1, 1], [1, 1]],
        'T': [[1, 1, 1], [0, 1, 0]],
        'skew': [[0, 1, 1], [1, 1, 0]],
        'L': [[1, 0], [1, 0], [1, 1]]
    }

    # The sequence of pieces that will enter the grid.
    PIECE_SEQUENCE = [
        "L", "T", "square", "T", "straight", "skew", "T", "square", "square", "skew", "straight", "square"
    ]

    # The target grid state we want to achieve.
    TARGET_GRID_STR = [
        "░░░░░░░░░░", "░░░░░░░░░░", "░░░░░░░░░░", "░░░░░░░░░░",
        "█░░░░░░░██", "███████░██", "░█████████", "██████████",
        "█░██░█████", "█████░████"
    ]
    TARGET_GRID = [[FILLED_CELL_VAL if char == '█' else EMPTY_CELL_VAL for char in row] for row in TARGET_GRID_STR]

    # The answer choices, each a long string of commands.
    ANSWERS = {
        "A": "shift left 4, down, shift left 1, rotate, rotate, down, shift right 2, down, down, rotate, shift right 5, down, rotate, shift right 2, down, shift left 3, down, down, shift right 1, down, rotate, shift left 4, down, rotate, shift right 4, down, shift right 4, down",
        "B": "shift left 4, down, shift left 1, rotate, rotate, down, shift right 2, down, down, rotate, shift right 5, down, shift right 2, down, shift left 3, down, down, shift left 2, down, rotate, shift left 4, down, rotate, shift right 4, down, shift right 4, down",
        "C": "shift left 4, down, shift left 1, rotate, rotate, down, shift right 2, down, down, rotate, shift right 5, down, rotate, shift right 2, down, shift left 2, down, down, shift left 2, down, rotate, shift left 4, down, rotate, shift right 4, down, shift right 4, down",
        "D": "shift left 4, down, shift left 1, rotate, rotate, down, shift right 4, down, down, rotate, shift right 5, down, rotate, shift right 2, down, shift left 3, down, down, shift left 2, down, rotate, shift left 4, down, rotate, shift right 4, down, shift right 4, down",
        "E": "shift left 4, down, shift left 2, rotate, rotate, down, shift right 2, down, down, rotate, shift right 5, down, rotate, shift right 2, down, shift left 3, down, down, shift left 2, down, rotate, shift left 4, down, rotate, shift right 4, down, shift right 4, down",
        "F": "shift left 4, down, shift left 1, rotate, rotate, down, shift right 2, down, down, rotate, shift right 5, down, rotate, shift right 2, down, shift left 3, down, down, shift left 2, down, rotate, shift left 4, down, rotate, shift right 4, down, shift right 4, down"
    }

    # --- Helper Functions ---
    def create_grid():
        """Creates a new empty grid."""
        return [[EMPTY_CELL_VAL for _ in range(GRID_WIDTH)] for _ in range(GRID_HEIGHT)]

    def rotate_piece(piece):
        """Rotates a piece matrix 90 degrees clockwise."""
        transposed = list(zip(*piece))
        return [list(row)[::-1] for row in transposed]

    def check_collision(grid, piece, pos):
        """Checks if a piece at a given position collides with walls or other blocks."""
        pos_r, pos_c = pos
        for r_offset, row_data in enumerate(piece):
            for c_offset, cell in enumerate(row_data):
                if cell == FILLED_CELL_VAL:
                    grid_r, grid_c = pos_r + r_offset, pos_c + c_offset
                    if not (0 <= grid_c < GRID_WIDTH and grid_r < GRID_HEIGHT):
                        return True  # Collision with walls or floor
                    if grid_r >= 0 and grid[grid_r][grid_c] == FILLED_CELL_VAL:
                        return True  # Collision with another piece
        return False

    def place_piece(grid, piece, pos):
        """Permanently places a piece onto the grid."""
        pos_r, pos_c = pos
        for r_offset, row_data in enumerate(piece):
            for c_offset, cell in enumerate(row_data):
                if cell == FILLED_CELL_VAL:
                    grid_r, grid_c = pos_r + r_offset, pos_c + c_offset
                    if 0 <= grid_r < GRID_HEIGHT and 0 <= grid_c < GRID_WIDTH:
                        grid[grid_r][grid_c] = FILLED_CELL_VAL

    def parse_commands(command_string):
        """Parses the command string into a list of command lists for each piece."""
        all_cmds = [c.strip() for c in command_string.split(',')]
        command_groups = []
        current_group = []
        for cmd in all_cmds:
            if cmd == 'down':
                command_groups.append(current_group)
                current_group = []
            else:
                current_group.append(cmd)
        return command_groups

    # --- Main Simulation Logic ---
    def run_simulation(command_string):
        """Runs a full game simulation for a given command string."""
        grid = create_grid()
        pieces = copy.deepcopy(PIECE_SEQUENCE)
        command_groups = parse_commands(command_string)

        if len(command_groups) != len(pieces):
            return None  # Should not happen with the given data

        for i in range(len(pieces)):
            piece_name = pieces[i]
            commands = command_groups[i]
            
            current_piece = copy.deepcopy(TETROMINOES[piece_name])
            
            # Start piece with its top-most solid part at row 0
            start_row_offset = 0
            for r_idx, row in enumerate(current_piece):
                if any(cell == 1 for cell in row):
                    start_row_offset = r_idx
                    break
            current_pos = [-start_row_offset, START_COL]

            # Apply pre-drop commands (rotate, shift)
            for cmd in commands:
                if cmd == "rotate":
                    current_piece = rotate_piece(current_piece)
                elif "shift left" in cmd:
                    amount = int(cmd.split()[-1])
                    current_pos[1] -= amount
                elif "shift right" in cmd:
                    amount = int(cmd.split()[-1])
                    current_pos[1] += amount
            
            # Drop the piece until it collides
            drop_pos_r = current_pos[0]
            while not check_collision(grid, current_piece, (drop_pos_r + 1, current_pos[1])):
                drop_pos_r += 1
            
            # Place the piece at its final position
            place_piece(grid, current_piece, (drop_pos_r, current_pos[1]))
        
        return grid

    # --- Find the Correct Answer ---
    for key, cmd_str in ANSWERS.items():
        final_grid = run_simulation(cmd_str)
        if final_grid == TARGET_GRID:
            print(key)
            return

solve_tetris_puzzle()