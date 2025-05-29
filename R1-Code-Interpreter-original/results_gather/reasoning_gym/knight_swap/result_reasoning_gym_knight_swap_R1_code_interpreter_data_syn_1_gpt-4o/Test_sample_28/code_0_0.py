from collections import deque

# Define the board and initial positions
board = {
    'A2': 'w', 'B1': 'w',
    'C3': 'B', 'D1': 'B',
    'A1': '.', 'B3': '.', 'C1': '.', 'D3': '.'
}

# Target positions
target_positions = {
    'w': ['C3', 'D1'],
    'B': ['A2', 'B1']
}

# Knight moves (L-shape)
knight_moves = [
    (2, 1), (2, -1), (-2, 1), (-2, -1),
    (1, 2), (1, -2), (-1, 2), (-1, -2)
]

# Convert board position to coordinates
def pos_to_coord(pos):
    col, row = pos
    return ord(col) - ord('A'), int(row) - 1

# Convert coordinates to board position
def coord_to_pos(coord):
    col, row = coord
    return chr(col + ord('A')) + str(row + 1)

# Check if a move is within bounds and to an empty square
def is_valid_move(start, end, board):
    if end not in board or board[end] != '.':
        return False
    start_coord = pos_to_coord(start)
    end_coord = pos_to_coord(end)
    dx, dy = end_coord[0] - start_coord[0], end_coord[1] - start_coord[1]
    return (dx, dy) in knight_moves

# BFS to find the sequence of moves
def find_moves():
    queue = deque([(board, [], 'w')])  # (current board, moves, current player)
    visited = set()

    while queue:
        current_board, moves, current_player = queue.popleft()
        board_tuple = tuple(sorted(current_board.items()))
        if board_tuple in visited:
            continue
        visited.add(board_tuple)

        # Check if the current board matches the target positions
        if all(current_board[pos] == 'w' for pos in target_positions['w']) and \
           all(current_board[pos] == 'B' for pos in target_positions['B']):
            return moves

        # Generate all possible moves for the current player
        for pos, piece in current_board.items():
            if piece == current_player:
                for dx, dy in knight_moves:
                    start_coord = pos_to_coord(pos)
                    end_coord = (start_coord[0] + dx, start_coord[1] + dy)
                    end_pos = coord_to_pos(end_coord)
                    if is_valid_move(pos, end_pos, current_board):
                        # Make the move
                        new_board = current_board.copy()
                        new_board[end_pos] = current_player
                        new_board[pos] = '.'
                        new_moves = moves + [f"{current_player},{pos},{end_pos}"]
                        # Alternate player
                        next_player = 'B' if current_player == 'w' else 'w'
                        queue.append((new_board, new_moves, next_player))

    return None

# Find the sequence of moves
result = find_moves()

# Output the result
if result:
    print(result)
else:
    print("No")