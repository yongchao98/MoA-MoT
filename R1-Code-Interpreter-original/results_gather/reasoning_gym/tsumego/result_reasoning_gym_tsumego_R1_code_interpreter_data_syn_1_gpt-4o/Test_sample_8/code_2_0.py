# Define the board
board = [
    ['.', '.', 'O', '.', '.', '.', '.', 'X', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'X', 'O', 'X', '.', '.', '.'],
    ['.', '.', '.', '.', '.', 'X', 'O', 'O', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'X', 'O', 'X', 'O', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', 'X', '.', '.', 'O', '.'],
    ['.', 'X', 'O', '.', 'O', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['O', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'X', '.', '.', '.'],
    ['.', '.', '.', 'X', '.', '.', '.', '.', 'O', '.', '.', '.'],
    ['.', '.', 'X', '.', '.', '.', '.', '.', '.', '.', '.', '.']
]

# Function to check if a move captures any stones
def is_capture_move(board, move, color):
    x, y = move
    opponent = 'O' if color == 'X' else 'X'
    board[y][x] = color  # Place the stone
    # Check if any opponent groups are captured
    def has_liberty(x, y, visited):
        if (x, y) in visited:
            return False
        visited.add((x, y))
        if board[y][x] == '.':
            return True
        if board[y][x] != opponent:
            return False
        # Check all adjacent positions
        return (x > 0 and has_liberty(x - 1, y, visited)) or \
               (x < len(board[0]) - 1 and has_liberty(x + 1, y, visited)) or \
               (y > 0 and has_liberty(x, y - 1, visited)) or \
               (y < len(board) - 1 and has_liberty(x, y + 1, visited))
    
    # Check all opponent stones
    for i in range(len(board)):
        for j in range(len(board[0])):
            if board[i][j] == opponent:
                if not has_liberty(j, i, set()):
                    return True
    return False

# Check if placing a stone at H2 captures any stones
move_h2 = (7, 10)  # H2
capture_h2 = is_capture_move(board, move_h2, 'X')
print(capture_h2)