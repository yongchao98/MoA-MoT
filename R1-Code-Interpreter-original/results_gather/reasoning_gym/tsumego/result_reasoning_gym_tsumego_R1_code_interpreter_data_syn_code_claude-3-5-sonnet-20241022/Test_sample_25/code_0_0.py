def is_valid_move(board, x, y):
    return 0 <= x < 9 and 0 <= y < 9 and board[y][x] == '.'

def get_liberties(board, x, y, visited):
    if not (0 <= x < 9 and 0 <= y < 9):
        return set()
    if (x, y) in visited:
        return set()
    if board[y][x] == '.':
        return {(x, y)}
    if board[y][x] != 'O':  # We're only checking white stones
        return set()
    
    visited.add((x, y))
    liberties = set()
    for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
        liberties.update(get_liberties(board, x + dx, y + dy, visited))
    return liberties

def can_capture_after_move(board, move_x, move_y):
    # Make a copy of the board and place the move
    new_board = [row[:] for row in board]
    new_board[move_y][move_x] = 'X'
    
    # Check all adjacent points for white stones
    for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
        x, y = move_x + dx, move_y + dy
        if 0 <= x < 9 and 0 <= y < 9 and new_board[y][x] == 'O':
            liberties = get_liberties(new_board, x, y, set())
            if not liberties:  # No liberties means capture is possible
                return True
    return False

# Initialize the board
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', 'O', '.'],
    ['.', '.', 'X', '.', '.', '.', '.', '.', 'O'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['X', 'X', 'X', 'X', '.', '.', 'X', 'X', '.'],
    ['.', 'X', 'O', 'O', 'X', 'O', 'X', '.', 'O'],
    ['.', 'X', 'O', '.', '.', '.', '.', '.', '.'],
    ['.', '.', 'O', 'X', '.', '.', '.', '.', '.']
]

# Check all empty points for capturing moves
for y in range(9):
    for x in range(9):
        if board[y][x] == '.' and can_capture_after_move(board, x, y):
            print(f"Capturing move found at column {chr(65+x)}, row {9-y}")