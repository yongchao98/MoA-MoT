from collections import deque

# Initial board setup
initial_board = {
    'A3': 'w', 'D2': 'w',
    'C2': 'B', 'C1': 'B',
    'B3': '.', 'B2': '.', 'A2': '.', 'B1': '.'
}

# Target board setup
target_board = {
    'A3': 'B', 'D2': 'B',
    'C2': 'w', 'C1': 'w',
    'B3': '.', 'B2': '.', 'A2': '.', 'B1': '.'
}

# Knight moves (L-shape)
knight_moves = [
    (2, 1), (2, -1), (-2, 1), (-2, -1),
    (1, 2), (1, -2), (-1, 2), (-1, -2)
]

# Convert position to coordinates
def pos_to_coord(pos):
    col, row = pos
    return ord(col) - ord('A'), int(row) - 1

# Convert coordinates to position
def coord_to_pos(x, y):
    return chr(x + ord('A')) + str(y + 1)

# Generate all possible moves for a knight
def generate_moves(pos):
    x, y = pos_to_coord(pos)
    moves = []
    for dx, dy in knight_moves:
        nx, ny = x + dx, y + dy
        if 0 <= nx < 4 and 0 <= ny < 3:
            new_pos = coord_to_pos(nx, ny)
            if new_pos in initial_board:  # Ensure the position is valid on the board
                moves.append(new_pos)
    return moves

# BFS to find the sequence of moves
def bfs_solve(initial, target):
    queue = deque([(initial, [])])
    visited = set()
    visited.add(tuple(sorted(initial.items())))

    while queue:
        current_board, moves = queue.popleft()

        if current_board == target:
            return moves

        for pos, knight in current_board.items():
            if knight in 'wB':
                for move in generate_moves(pos):
                    if current_board[move] == '.':
                        # Make the move
                        new_board = current_board.copy()
                        new_board[move] = knight
                        new_board[pos] = '.'
                        new_moves = moves + [f"{knight},{pos},{move}"]

                        board_tuple = tuple(sorted(new_board.items()))
                        if board_tuple not in visited:
                            visited.add(board_tuple)
                            queue.append((new_board, new_moves))

    return "No"

# Solve the puzzle
result = bfs_solve(initial_board, target_board)
print(result)