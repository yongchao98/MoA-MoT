# Define the board
board = [
    ['.', 'X', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['O', '.', '.', 'O', '.', '.', '.', '.', '.'],
    ['X', '.', '.', '.', '.', '.', '.', 'X', '.'],
    ['X', '.', '.', '.', '.', '.', 'X', 'O', '.'],
    ['.', '.', '.', 'X', '.', 'X', '.', '.', '.'],
    ['.', '.', '.', '.', 'O', 'O', 'X', '.', '.'],
    ['.', '.', '.', 'X', 'O', 'O', '.', 'X', '.'],
    ['.', '.', '.', '.', 'X', 'O', 'X', '.', '.']
]

# Function to check if a move captures stones
def check_capture(board, move):
    x, y = move
    board[y][x] = 'X'  # Place the black stone
    # Check the liberties of the white group
    white_group = [(3, 2), (4, 2), (3, 1), (4, 1)]
    for wx, wy in white_group:
        if board[wy][wx] == 'O':
            # Check liberties
            liberties = 0
            for nx, ny in [(wx-1, wy), (wx+1, wy), (wx, wy-1), (wx, wy+1)]:
                if 0 <= nx < 9 and 0 <= ny < 9 and board[ny][nx] == '.':
                    liberties += 1
            if liberties == 0:
                return True
    return False

# Check if playing at F3 captures the white stones
move = (5, 2)  # F3
capture = check_capture(board, move)
print(capture)