from collections import deque

# Define the board and initial positions
board = {
    'A3': 'B', 'B3': None, 'C3': '.', 'D3': 'B',
    'A2': '.', 'B2': None, 'C2': 'w', 'D2': 'w',
    'A1': '.', 'B1': '.', 'C1': '.', 'D1': None
}

# Define knight moves
knight_moves = [(2, 1), (2, -1), (-2, 1), (-2, -1), (1, 2), (1, -2), (-1, 2), (-1, -2)]

# Convert board position to coordinates
def pos_to_coord(pos):
    col, row = pos
    return ord(col) - ord('A'), int(row) - 1

# Convert coordinates to board position
def coord_to_pos(x, y):
    return chr(x + ord('A')) + str(y + 1)

# Check if a move is within bounds
def is_within_bounds(x, y):
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

        # Check if the goal is achieved
        if current_board['A3'] == 'w' and current_board['D3'] == 'w' and current_board['C2'] == 'B' and current_board['D2'] == 'B':
            return moves

        # Find all knights of the current player
        for pos, knight in current_board.items():
            if knight == current_player:
                x, y = pos_to_coord(pos)
                # Try all possible knight moves
                for dx, dy in knight_moves:
                    nx, ny = x + dx, y + dy
                    new_pos = coord_to_pos(nx, ny)
                    if is_within_bounds(nx, ny) and current_board.get(new_pos) == '.':
                        # Make the move
                        new_board = current_board.copy()
                        new_board[new_pos] = current_player
                        new_board[pos] = '.'
                        new_moves = moves + [f"{current_player},{pos},{new_pos}"]
                        # Alternate player
                        next_player = 'B' if current_player == 'w' else 'w'
                        queue.append((new_board, new_moves, next_player))

    return "No"

# Run the BFS to find the solution
result = bfs()
print(result)