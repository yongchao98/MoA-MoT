from collections import deque

# Define the board and initial positions
board = {
    'A1': None, 'B1': None, 'C1': 'B', 'D1': 'w',
    'A2': None, 'B2': 'w', 'C2': None, 'D2': None,
    'A3': None, 'B3': None, 'C3': 'B', 'D3': None
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
def coord_to_pos(coord):
    col, row = coord
    return chr(col + ord('A')) + str(row + 1)

# Check if a move is valid
def is_valid_move(start, end):
    if end not in board or board[end] is not None:
        return False
    start_x, start_y = pos_to_coord(start)
    end_x, end_y = pos_to_coord(end)
    for dx, dy in knight_moves:
        if start_x + dx == end_x and start_y + dy == end_y:
            return True
    return False

# BFS to find the sequence of moves
def find_moves():
    queue = deque([([], board)])
    visited = set()

    while queue:
        moves, current_board = queue.popleft()
        board_state = tuple(sorted(current_board.items()))
        if board_state in visited:
            continue
        visited.add(board_state)

        # Check if the current board is the target state
        if current_board['B2'] == 'B' and current_board['D1'] == 'B' and \
           current_board['C3'] == 'w' and current_board['C1'] == 'w':
            return moves

        # Generate possible moves
        for pos, knight in current_board.items():
            if knight is None:
                continue
            for dx, dy in knight_moves:
                new_x, new_y = pos_to_coord(pos)
                new_x += dx
                new_y += dy
                new_pos = coord_to_pos((new_x, new_y))
                if is_valid_move(pos, new_pos):
                    new_board = current_board.copy()
                    new_board[new_pos] = knight
                    new_board[pos] = None
                    new_moves = moves + [f"{knight},{pos},{new_pos}"]
                    queue.append((new_moves, new_board))

    return "No"

# Find and print the moves
result = find_moves()
print(result)