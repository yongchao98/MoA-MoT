# Define the board
board = [
    ['.', '.', '.', '.', '.', 'O', '.', '.', 'X', 'O', 'X', '.'],
    ['.', '.', '.', '.', '.', '.', '.', 'X', 'O', 'O', '.', 'X'],
    ['.', 'O', '.', '.', '.', 'X', 'X', '.', 'X', 'O', 'X', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'X', 'O', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'X'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', 'X'],
    ['.', '.', '.', '.', '.', '.', 'O', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', 'O', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', 'X', '.', '.', 'O', '.', '.']
]

# Check if placing a black stone at C2 captures the white stone at C3
def is_capturable(x, y):
    # Check if the position is within bounds
    if x < 0 or x >= len(board[0]) or y < 0 or y >= len(board):
        return False
    # Check if the position is a white stone
    if board[y][x] != 'O':
        return False
    # Check the liberties of the white stone
    liberties = 0
    for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
        nx, ny = x + dx, y + dy
        if 0 <= nx < len(board[0]) and 0 <= ny < len(board):
            if board[ny][nx] == '.':
                liberties += 1
    return liberties == 1

# Check if placing a black stone at C2 captures the white stone at C3
capture_move = (2, 1)  # C2
x, y = capture_move
if is_capturable(x, y + 1):  # Check the stone at C3
    print("<<<C2>>>")
else:
    print("No capture possible")