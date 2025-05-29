# Define the board
board = [
    ['.', '.', 'O', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'O', '.'],
    ['.', '.', 'X', '.', '.', '.', 'O', '.', 'X', '.'],
    ['.', '.', 'X', '.', '.', '.', 'X', 'X', '.', '.'],
    ['.', '.', '.', '.', 'O', 'X', 'O', 'O', 'X', '.'],
    ['.', '.', '.', '.', '.', 'X', 'O', '.', '.', '.'],
    ['.', 'X', '.', '.', '.', '.', 'X', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', 'O', '.', '.'],
    ['.', '.', 'X', '.', '.', '.', '.', '.', '.', '.']
]

# Function to check if a move at (x, y) captures any stones
def is_capture_move(x, y, board):
    # Check if the move is within bounds and the position is empty
    if 0 <= x < 10 and 0 <= y < 10 and board[y][x] == '.':
        # Temporarily place a black stone at (x, y)
        board[y][x] = 'X'
        
        # Check the group at F6, G6, F5
        group = [(5, 6), (6, 6), (5, 5)]
        for gx, gy in group:
            if board[gy][gx] == 'O':
                # Check liberties
                liberties = 0
                for nx, ny in [(gx-1, gy), (gx+1, gy), (gx, gy-1), (gx, gy+1)]:
                    if 0 <= nx < 10 and 0 <= ny < 10 and board[ny][nx] == '.':
                        liberties += 1
                if liberties == 0:
                    # Restore the board
                    board[y][x] = '.'
                    return True
        # Restore the board
        board[y][x] = '.'
    return False

# Check if E6 is a capture move
capture_move = is_capture_move(4, 5, board)
print(capture_move)