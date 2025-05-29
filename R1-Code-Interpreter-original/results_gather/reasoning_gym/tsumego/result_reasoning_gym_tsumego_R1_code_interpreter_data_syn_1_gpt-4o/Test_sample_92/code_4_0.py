# Define the board
board = [
    ['.', '.', '.', '.', '.', '.', 'X', '.', 'X', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'X', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'X', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['X', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.'],
    ['X', '.', 'X', '.', 'O', '.', '.', '.', '.', '.', '.'],
    ['.', '.', '.', '.', '.', '.', '.', '.', 'X', '.', '.'],
    ['X', '.', '.', '.', '.', '.', 'X', 'O', 'O', '.', '.'],
    ['.', '.', '.', '.', '.', '.', 'X', 'O', 'O', '.', 'X'],
    ['.', '.', '.', '.', 'O', '.', 'X', 'O', 'X', '.', '.']
]

# Manually check the move G2
def manual_capture_check(board):
    # Place the stone at G2
    board[9][6] = 'X'
    
    # Check if the group at G3, H3, H2 is surrounded
    captured = 0
    if board[8][7] == 'O' and board[9][7] == 'O' and board[9][8] == 'O':
        if board[8][6] == 'X' and board[9][5] == 'X' and board[10][7] == 'X' and board[9][9] == 'X':
            captured = 3  # G3, H3, H2 are captured
    
    return captured

captured_stones = manual_capture_check(board)
print(captured_stones)