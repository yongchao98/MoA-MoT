# Define the initial and target positions
initial_positions = {'w': ['A2', 'B1'], 'B': ['C3', 'C2']}
target_positions = {'w': ['C3', 'C2'], 'B': ['A2', 'B1']}

# Define the board and knight moves
board = [['.', '.', 'B', '.'],
         ['.', '.', 'B', 'w'],
         ['.', 'w', '.', '.']]

# Convert board positions to coordinates
def pos_to_coord(pos):
    col, row = pos
    return ord(col) - ord('A'), int(row) - 1

# Convert coordinates to board positions
def coord_to_pos(coord):
    col, row = coord
    return chr(col + ord('A')) + str(row + 1)

# Possible knight moves
knight_moves = [(2, 1), (2, -1), (-2, 1), (-2, -1),
                (1, 2), (1, -2), (-1, 2), (-1, -2)]

# Check if a move is valid
def is_valid_move(start, end, board):
    x, y = start
    nx, ny = end
    if 0 <= nx < 4 and 0 <= ny < 3 and board[ny][nx] == '.':
        return True
    return False

# Simulate the moves
def simulate_moves():
    moves = []
    # Start with black's move
    current_color = 'B'
    while initial_positions['w'] != target_positions['w'] or initial_positions['B'] != target_positions['B']:
        moved = False
        for i, pos in enumerate(initial_positions[current_color]):
            x, y = pos_to_coord(pos)
            for dx, dy in knight_moves:
                nx, ny = x + dx, y + dy
                if is_valid_move((x, y), (nx, ny), board):
                    # Make the move
                    board[y][x] = '.'
                    board[ny][nx] = current_color
                    new_pos = coord_to_pos((nx, ny))
                    initial_positions[current_color][i] = new_pos
                    moves.append(f"{current_color},{pos},{new_pos}")
                    moved = True
                    break
            if moved:
                break
        if not moved:
            return "No"
        # Alternate turns
        current_color = 'w' if current_color == 'B' else 'B'
    return moves

# Get the result
result = simulate_moves()
print(result)