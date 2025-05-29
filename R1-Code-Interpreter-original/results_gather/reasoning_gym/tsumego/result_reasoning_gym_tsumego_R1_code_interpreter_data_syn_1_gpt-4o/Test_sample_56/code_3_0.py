# Define the board
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'X'],
    ['.', 'X', '.', '.', '.', '.', '.', 'X', '.', '.', '.', '.'],
    ['.', '.', 'O', '.', '.', 'O', '.', 'O', '.', 'X', 'X', '.'],
    ['.', '.', '.', '.', '.', '.', '.', 'X', 'X', 'O', 'O', 'X'],
    ['X', '.', '.', '.', '.', '.', '.', '.', 'X', 'O', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['O', '.', '.', '.', '.', '.', 'X', 'X', '.', '.', '.', '.'],
    ['O', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.']
]

# Check the liberties of a stone at a given position
def check_liberties(x, y):
    liberties = 0
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if 0 <= nx < 12 and 0 <= ny < 12 and board[nx][ny] == '.':
            liberties += 1
    return liberties

# Check if playing at G8 captures the stones at H8 and I8
def is_capture_move():
    # Simulate the move at G8 (row 4, column 6)
    board[4][6] = 'X'
    # Check liberties of the group at H8 and I8
    group_liberties = check_liberties(4, 7) + check_liberties(4, 8)
    # Restore the board
    board[4][6] = '.'
    return group_liberties == 0

# Check if G8 is a capturing move
capture_move = is_capture_move()
print(capture_move)