# Define the board
board = [
    ['.', '.', '.', '.', '.', '.', '.', '.', 'O', '.', '.', '.'],
    ['X', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', 'X', '.', 'X', 'X', '.', '.', '.', '.'],
    ['.', '.', '.', '.', 'O', 'X', 'O', 'O', 'X', '.', '.', 'X'],
    ['.', 'X', '.', '.', '.', 'X', 'O', '.', '.', '.', '.', '.'],
    ['.', 'X', 'X', '.', '.', '.', 'X', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['O', '.', '.', '.', '.', 'X', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'O', '.', '.', 'X', '.', '.']
]

# Function to check liberties of a stone
def check_liberties(x, y, color):
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    liberties = 0
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if 0 <= nx < 12 and 0 <= ny < 12 and board[nx][ny] == '.':
            liberties += 1
    return liberties

# Function to find the best move for Black
def find_best_move():
    max_captures = 0
    best_move = None
    for x in range(12):
        for y in range(12):
            if board[x][y] == '.':
                # Simulate placing a Black stone
                board[x][y] = 'X'
                captures = 0
                # Check surrounding White stones
                for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
                    nx, ny = x + dx, y + dy
                    if 0 <= nx < 12 and 0 <= ny < 12 and board[nx][ny] == 'O':
                        if check_liberties(nx, ny, 'O') == 0:
                            captures += 1
                # Update best move if more captures
                if captures > max_captures:
                    max_captures = captures
                    best_move = (x, y)
                # Revert the move
                board[x][y] = '.'
    return best_move

# Find the best move
best_move = find_best_move()
# Convert to board coordinates
if best_move:
    row, col = best_move
    col_letter = chr(ord('A') + col)
    row_number = 12 - row
    print(f'<<<{col_letter}{row_number}>>>')
else:
    print("No capturing move found")