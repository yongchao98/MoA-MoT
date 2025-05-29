from collections import deque

# Define the board and initial positions
board = {
    'A3': 'w', 'B1': 'w',
    'C3': 'B', 'D2': 'B',
    'B3': '.', 'D1': '.'
}

# Define the target positions
target_positions = {
    'w': ['C3', 'D2'],
    'B': ['A3', 'B1']
}

# Define knight moves
knight_moves = [
    (2, 1), (2, -1), (-2, 1), (-2, -1),
    (1, 2), (1, -2), (-1, 2), (-1, -2)
]

# Convert board positions to coordinates
def pos_to_coord(pos):
    col, row = pos
    return ord(col) - ord('A'), int(row) - 1

# Convert coordinates to board positions
def coord_to_pos(coord):
    col, row = coord
    return chr(col + ord('A')) + str(row + 1)

# Check if a position is within the board
def is_within_board(coord):
    col, row = coord
    return 0 <= col < 4 and 0 <= row < 3

# BFS to find the sequence of moves
def bfs_knight_swap():
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

        # Get current player's knight positions
        knight_positions = [pos for pos, piece in current_board.items() if piece == current_player]

        # Try all possible moves for each knight
        for pos in knight_positions:
            col, row = pos_to_coord(pos)
            for move in knight_moves:
                new_col, new_row = col + move[0], row + move[1]
                new_pos = coord_to_pos((new_col, new_row))
                
                if is_within_board((new_col, new_row)) and current_board.get(new_pos, '.') == '.':
                    # Make the move
                    new_board = current_board.copy()
                    new_board[new_pos] = current_player
                    new_board[pos] = '.'
                    new_moves = moves + [f"{current_player},{pos},{new_pos}"]
                    next_player = 'B' if current_player == 'w' else 'w'
                    queue.append((new_board, new_moves, next_player))

    return "No"

# Find the solution
solution = bfs_knight_swap()
print(solution)