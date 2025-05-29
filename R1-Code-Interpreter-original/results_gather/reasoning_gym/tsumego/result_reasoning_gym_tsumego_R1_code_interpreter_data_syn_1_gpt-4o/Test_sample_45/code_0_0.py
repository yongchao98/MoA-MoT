# Define the board
board = [
    ['.', '.', '.', '.', 'X', '.', '.', '.', '.', '.', '.'],
    ['.', 'O', '.', '.', 'X', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', '.'],
    ['.', '.', '.', 'X', '.', '.', '.', '.', 'O', 'X', 'X'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'X', 'O', 'O'],
    ['.', '.', '.', 'X', '.', '.', '.', '.', 'X', 'O', '.'],
    ['.', '.', 'X', '.', '.', '.', '.', '.', '.', 'X', '.'],
    ['.', '.', 'O', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', 'X', 'O'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.']
]

# Function to check if a position is within the board
def is_within_board(x, y):
    return 0 <= x < 11 and 0 <= y < 11

# Function to get the liberties of a stone
def get_liberties(x, y, color):
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    liberties = 0
    for dx, dy in directions:
        nx, ny = x + dx, y + dy
        if is_within_board(nx, ny) and board[nx][ny] == '.':
            liberties += 1
    return liberties

# Find the best move for Black
best_move = None
max_captures = 0

for x in range(11):
    for y in range(11):
        if board[x][y] == 'O':
            liberties = get_liberties(x, y, 'O')
            if liberties == 1:
                # Find the liberty position
                for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
                    nx, ny = x + dx, y + dy
                    if is_within_board(nx, ny) and board[nx][ny] == '.':
                        # Check if placing a black stone here captures the white stone
                        captures = 1  # Capturing this stone
                        if captures > max_captures:
                            max_captures = captures
                            best_move = (nx, ny)

# Convert the best move to board coordinates
if best_move:
    row, col = best_move
    col_letter = chr(ord('A') + col)
    row_number = 11 - row
    print(f'<<<{col_letter}{row_number}>>>')
else:
    print("No immediate capture move found.")