from collections import deque

# Define the board and initial positions
board = {
    'A1': None, 'B1': 'w', 'C1': None, 'D1': None,
    'A2': None, 'B2': None, 'C2': 'w', 'D2': 'B',
    'A3': None, 'B3': None, 'C3': 'B', 'D3': None
}

# Define the target positions
target_positions = {
    'w': ['C3', 'D2'],
    'B': ['C2', 'B1']
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
def bfs():
    queue = deque()
    queue.append((board, [], 'w'))  # (current board, moves, current player)
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

        # Get positions of current player's knights
        knight_positions = [pos for pos, knight in current_board.items() if knight == current_player]

        for pos in knight_positions:
            x, y = pos_to_coord(pos)
            for dx, dy in knight_moves:
                nx, ny = x + dx, y + dy
                if is_within_board(nx, ny):
                    new_pos = coord_to_pos(nx, ny)
                    if current_board[new_pos] is None:  # Move to empty square
                        new_board = current_board.copy()
                        new_board[new_pos] = current_player
                        new_board[pos] = None
                        new_moves = moves + [f"{current_player},{pos},{new_pos}"]
                        next_player = 'B' if current_player == 'w' else 'w'
                        queue.append((new_board, new_moves, next_player))

    return None

# Find the solution
solution = bfs()
if solution:
    print(solution)
else:
    print("No")