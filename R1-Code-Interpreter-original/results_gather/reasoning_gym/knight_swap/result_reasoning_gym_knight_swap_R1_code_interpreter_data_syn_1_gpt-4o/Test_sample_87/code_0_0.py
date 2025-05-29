from collections import deque

# Define the board and initial positions
board = {
    'A3': 'B', 'B3': 'w', 'C3': '.', 'D3': '.',
    'A2': 'w', 'B2': '.', 'C2': '.', 'D2': 'B',
    'A1': '.', 'B1': '.', 'C1': '.', 'D1': '.'
}

# Define target positions
target_positions = {
    'w': ['A3', 'D2'],
    'B': ['B3', 'A2']
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

# Generate possible moves for a knight
def generate_moves(pos):
    x, y = pos_to_coord(pos)
    moves = []
    for dx, dy in knight_moves:
        nx, ny = x + dx, y + dy
        if is_within_board(nx, ny):
            new_pos = coord_to_pos(nx, ny)
            if board[new_pos] == '.':
                moves.append(new_pos)
    return moves

# BFS to find the sequence of moves
def bfs():
    queue = deque([(board, [], 'w')])  # (current board, moves, next turn)
    visited = set()

    while queue:
        current_board, moves, turn = queue.popleft()
        board_tuple = tuple(sorted(current_board.items()))
        if board_tuple in visited:
            continue
        visited.add(board_tuple)

        # Check if the current board matches the target positions
        if all(current_board[pos] == 'w' for pos in target_positions['w']) and \
           all(current_board[pos] == 'B' for pos in target_positions['B']):
            return moves

        # Generate moves for the current turn
        for pos, knight in current_board.items():
            if knight == turn:
                for new_pos in generate_moves(pos):
                    new_board = current_board.copy()
                    new_board[new_pos] = knight
                    new_board[pos] = '.'
                    new_moves = moves + [f"{knight},{pos},{new_pos}"]
                    next_turn = 'B' if turn == 'w' else 'w'
                    queue.append((new_board, new_moves, next_turn))

    return None

# Find the solution
solution = bfs()
if solution:
    print(solution)
else:
    print("No")