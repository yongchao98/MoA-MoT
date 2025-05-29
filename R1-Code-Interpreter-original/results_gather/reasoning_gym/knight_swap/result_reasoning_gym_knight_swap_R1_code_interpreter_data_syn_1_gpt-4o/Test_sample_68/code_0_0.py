from collections import deque

# Define the board and initial positions
board = {
    'A1': 'B', 'B1': '.', 'C1': 'w', 'D1': '.',
    'A2': '.', 'B2': '.', 'C2': '.', 'D2': 'B',
    'A3': '.', 'B3': '.', 'C3': 'w', 'D3': '.'
}

# Define the target positions
target_positions = {
    'w': ['D2', 'A1'],
    'B': ['C3', 'C1']
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
    queue = deque([(board, [], 'w')])  # (current board, moves, current player)
    visited = set()

    while queue:
        current_board, moves, current_player = queue.popleft()
        board_tuple = tuple(current_board.items())

        if board_tuple in visited:
            continue
        visited.add(board_tuple)

        # Check if the current board matches the target positions
        if all(current_board[pos] == 'w' for pos in target_positions['w']) and \
           all(current_board[pos] == 'B' for pos in target_positions['B']):
            return moves

        # Get current player's knights
        for pos, knight in current_board.items():
            if knight == current_player:
                x, y = pos_to_coord(pos)

                # Try all possible knight moves
                for dx, dy in knight_moves:
                    nx, ny = x + dx, y + dy
                    new_pos = coord_to_pos(nx, ny)

                    if is_within_board(nx, ny) and current_board[new_pos] == '.':
                        # Make the move
                        new_board = current_board.copy()
                        new_board[new_pos] = current_player
                        new_board[pos] = '.'

                        # Add the move to the list
                        new_moves = moves + [f"{current_player},{pos},{new_pos}"]

                        # Alternate player
                        next_player = 'B' if current_player == 'w' else 'w'

                        # Add to queue
                        queue.append((new_board, new_moves, next_player))

    return "No"

# Find the solution
solution = bfs()
print(solution)