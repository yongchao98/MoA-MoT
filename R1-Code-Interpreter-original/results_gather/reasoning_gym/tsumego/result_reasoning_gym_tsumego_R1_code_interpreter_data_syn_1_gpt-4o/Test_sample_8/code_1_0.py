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

# Check if placing a stone at C6 or D6 captures any stones
move_c6 = (2, 6)  # C6
move_d6 = (3, 6)  # D6
capture_c6 = is_capture_move(board, move_c6, 'X')
capture_d6 = is_capture_move(board, move_d6, 'X')
print(capture_c6, capture_d6)