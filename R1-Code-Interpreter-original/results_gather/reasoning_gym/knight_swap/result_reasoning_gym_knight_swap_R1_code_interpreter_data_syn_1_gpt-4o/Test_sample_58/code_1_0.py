from collections import deque

# Define the initial board state with all positions
initial_board = {
    'A3': 'B', 'B3': 'w', 'C3': '.', 'D3': 'B',
    'A2': 'w', 'B2': '.', 'C2': '.', 'D2': '.',
    'A1': '.', 'B1': '.', 'C1': '.', 'D1': '.'
}

# Define the target board state
target_board = {
    'A3': 'w', 'B3': 'B', 'C3': '.', 'D3': 'w',
    'A2': 'B', 'B2': '.', 'C2': '.', 'D2': '.',
    'A1': '.', 'B1': '.', 'C1': '.', 'D1': '.'
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

# BFS to find the sequence of moves
def bfs_knight_swap(initial, target):
    queue = deque([(initial, [], 'B')])  # (board_state, move_sequence, current_turn)
    visited = set()

    while queue:
        board, moves, turn = queue.popleft()
        board_tuple = tuple(sorted(board.items()))
        
        if board_tuple in visited:
            continue
        visited.add(board_tuple)

        if board == target:
            return moves

        # Get all knights of the current turn
        knights = [pos for pos, color in board.items() if color == turn]

        for knight in knights:
            x, y = pos_to_coord(knight)
            for dx, dy in knight_moves:
                nx, ny = x + dx, y + dy
                if is_within_board(nx, ny):
                    new_pos = coord_to_pos(nx, ny)
                    if board[new_pos] == '.':
                        # Make the move
                        new_board = board.copy()
                        new_board[new_pos] = turn
                        new_board[knight] = '.'
                        new_moves = moves + [f"{turn},{knight},{new_pos}"]
                        # Alternate turn
                        next_turn = 'w' if turn == 'B' else 'B'
                        queue.append((new_board, new_moves, next_turn))

    return "No"

# Find the sequence of moves
result = bfs_knight_swap(initial_board, target_board)
print(result)