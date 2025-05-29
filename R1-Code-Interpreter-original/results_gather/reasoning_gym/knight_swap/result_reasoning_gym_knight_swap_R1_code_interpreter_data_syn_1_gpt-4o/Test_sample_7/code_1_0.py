# Define the board and initial positions
board = {
    'A1': None, 'B1': None, 'C1': 'B', 'D1': None,
    'A2': None, 'B2': None, 'C2': 'w', 'D2': 'w',
    'A3': None, 'B3': 'B', 'C3': None, 'D3': None
}

# Define the target positions for each knight
target_positions = {
    'w': ['B3', 'C1'],
    'B': ['C2', 'D2']
}

# Define knight moves
knight_moves = [
    (2, 1), (2, -1), (-2, 1), (-2, -1),
    (1, 2), (1, -2), (-1, 2), (-1, -2)
]

# Convert board position to coordinates
def pos_to_coord(pos):
    col, row = pos
    return ord(col) - ord('A'), int(row) - 1

# Convert coordinates to board position
def coord_to_pos(x, y):
    return chr(x + ord('A')) + str(y + 1)

# Check if a position is within the board
def is_within_board(x, y):
    return 0 <= x < 4 and 0 <= y < 3

# DFS to find the sequence of moves with a depth limit
def dfs(board, moves, current_player, depth_limit):
    if depth_limit == 0:
        return None

    # Check if current board is a solution
    if all(board[pos] == 'w' for pos in target_positions['w']) and \
       all(board[pos] == 'B' for pos in target_positions['B']):
        return moves

    current_positions = [pos for pos, piece in board.items() if piece == current_player]

    # Generate possible moves
    for pos in current_positions:
        x, y = pos_to_coord(pos)
        for dx, dy in knight_moves:
            nx, ny = x + dx, y + dy
            new_pos = coord_to_pos(nx, ny)
            if is_within_board(nx, ny) and board[new_pos] is None:
                # Make the move
                new_board = board.copy()
                new_board[new_pos] = current_player
                new_board[pos] = None
                new_moves = moves + [f"{current_player},{pos},{new_pos}"]
                result = dfs(new_board, new_moves, 'B' if current_player == 'w' else 'w', depth_limit - 1)
                if result is not None:
                    return result

    return None

# Try to find a solution with increasing depth limits
for depth_limit in range(1, 15):
    solution = dfs(board, [], 'w', depth_limit)
    if solution:
        print(solution)
        break
else:
    print("No")