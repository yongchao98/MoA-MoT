def solve():
    """
    This function solves the Tetris-like puzzle by simulating the game for each of the provided command options.
    It identifies which sequence of commands results in the given target grid configuration.
    """

    # --- Game Definitions ---

    PIECES = {
        'straight': [[1, 1, 1, 1]],
        'square': [[1, 1], [1, 1]],
        'T': [[1, 1, 1], [0, 1, 0]],
        'skew': [[0, 1, 1], [1, 1, 0]],
        'L': [[1, 0], [1, 0], [1, 1]]
    }

    PIECE_SEQUENCE = [
        "L", "T", "square", "T", "straight", "skew", "T", "square", "square", "skew", "straight", "square"
    ]

    TARGET_GRID = [
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 0, 0, 0, 0, 0, 0, 1, 1],
        [1, 1, 1, 1, 1, 1, 1, 0, 1, 1],
        [0, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        [1, 0, 1, 1, 0, 1, 1, 1, 1, 1],
        [1, 1, 1, 1, 1, 0, 1, 1, 1, 1]
    ]

    COMMANDS = {
        'A': "shift left 4, down, shift left 1, rotate, rotate, down, shift right 2, down, down, rotate, shift right 5, down, rotate, shift right 2, down, shift left 3, down, down, shift right 1, down, rotate, shift left 4, down, rotate, shift right 4, down, shift right 4, down",
        'B': "shift left 4, down, shift left 1, rotate, rotate, down, shift right 2, down, down, rotate, shift right 5, down, shift right 2, down, shift left 3, down, down, shift left 2, down, rotate, shift left 4, down, rotate, shift right 4, down, shift right 4, down",
        'C': "shift left 4, down, shift left 1, rotate, rotate, down, shift right 2, down, down, rotate, shift right 5, down, rotate, shift right 2, down, shift left 2, down, down, shift left 2, down, rotate, shift left 4, down, rotate, shift right 4, down, shift right 4, down",
        'D': "shift left 4, down, shift left 1, rotate, rotate, down, shift right 4, down, down, rotate, shift right 5, down, rotate, shift right 2, down, shift left 3, down, down, shift left 2, down, rotate, shift left 4, down, rotate, shift right 4, down, shift right 4, down",
        'E': "shift left 4, down, shift left 2, rotate, rotate, down, shift right 2, down, down, rotate, shift right 5, down, rotate, shift right 2, down, shift left 3, down, down, shift left 2, down, rotate, shift left 4, down, rotate, shift right 4, down, shift right 4, down",
        'F': "shift left 4, down, shift left 1, rotate, rotate, down, shift right 2, down, down, rotate, shift right 5, down, rotate, shift right 2, down, shift left 3, down, down, shift left 2, down, rotate, shift left 4, down, rotate, shift right 4, down, shift right 4, down"
    }

    class TetrisSimulator:
        def __init__(self, pieces_map):
            self.GRID_WIDTH = 10
            self.GRID_HEIGHT = 10
            self.SPAWN_X = 4
            self.PIECES = pieces_map
            self.grid = []
            self.current_piece = {}

        def _create_empty_grid(self):
            return [[0 for _ in range(self.GRID_WIDTH)] for _ in range(self.GRID_HEIGHT)]

        def _rotate(self, shape):
            return [list(row) for row in zip(*shape[::-1])]

        def _check_collision(self, piece_shape, piece_x, piece_y):
            for r_idx, row in enumerate(piece_shape):
                for c_idx, cell in enumerate(row):
                    if cell:
                        grid_x = piece_x + c_idx
                        grid_y = piece_y + r_idx
                        if not (0 <= grid_x < self.GRID_WIDTH and 0 <= grid_y < self.GRID_HEIGHT):
                            return True
                        if self.grid[grid_y][grid_x] == 1:
                            return True
            return False

        def _place_piece(self):
            for r_idx, row in enumerate(self.current_piece['shape']):
                for c_idx, cell in enumerate(row):
                    if cell:
                        self.grid[self.current_piece['y'] + r_idx][self.current_piece['x'] + c_idx] = 1
        
        def _parse_command_string(self, command_string):
            commands_per_piece = []
            current_piece_commands = []
            for cmd_part in command_string.split(', '):
                parts = cmd_part.split()
                current_piece_commands.append(parts)
                if parts[0] == 'down':
                    commands_per_piece.append(current_piece_commands)
                    current_piece_commands = []
            return commands_per_piece

        def run(self, piece_sequence, command_string):
            self.grid = self._create_empty_grid()
            all_commands = self._parse_command_string(command_string)

            for piece_name, commands in zip(piece_sequence, all_commands):
                self.current_piece = {'shape': [list(row) for row in self.PIECES[piece_name]], 'x': self.SPAWN_X, 'y': 0}

                for cmd in commands:
                    cmd_type = cmd[0]
                    if cmd_type == 'shift':
                        direction = cmd[1]
                        amount = int(cmd[2])
                        self.current_piece['x'] += -amount if direction == 'left' else amount
                    elif cmd_type == 'rotate':
                        self.current_piece['shape'] = self._rotate(self.current_piece['shape'])
                    elif cmd_type == 'down':
                        # The piece starts at y=0, check if it's a valid starting position before dropping
                        if self._check_collision(self.current_piece['shape'], self.current_piece['x'], self.current_piece['y']):
                             # Invalid initial placement, something is wrong
                             return None

                        while True:
                            if self._check_collision(self.current_piece['shape'], self.current_piece['x'], self.current_piece['y'] + 1):
                                self._place_piece()
                                break
                            else:
                                self.current_piece['y'] += 1
            return self.grid

    simulator = TetrisSimulator(PIECES)
    final_answer = "No solution found"
    for option, command_string in COMMANDS.items():
        final_grid = simulator.run(PIECE_SEQUENCE, command_string)
        if final_grid == TARGET_GRID:
            final_answer = option
            break
            
    print(f'<<<{final_answer}>>>')

solve()